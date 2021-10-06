using BenchmarkTools
using JuMP, NLPModelsJuMP, NLPModelsAlgencan
import LazySets: Ellipsoid
import Base: in
"""
Structure of an Ellipsoid satisfying dot(x,A*x) + 2*dot(b,x) ≤ α
"""
@kwdef struct EllipsoidBBSI
    A::SparseMatrixCSC
    b::Vector
    α::Float64
end


function in(x₀::Vector, ell::EllipsoidBBSI)
    A, b, α  = ell.A, ell.b, ell.α
    if dot(x₀,A*x₀) + 2*dot(b,x₀) ≤ α
        return true
    else
        return false
    end
end
"""
Transform Ellipsoid in format dot(x-c,Q⁻¹*(x-c)) ≤ 1 
into format  dot(x,A*x) + 2*dot(b,x) ≤ α
from shape matrix Q and center of ellipsoid c
"""
function EllipsoidBBSI(c::Vector, Q::AbstractMatrix)
    A = inv(Matrix(Q))
    b = -A*c
    α = 1 + dot(c,b)
    return EllipsoidBBSI(A,b,α)
end


"""
Transform Ellipsoid in format dot(x-c,Q⁻¹*(x-c)) ≤ 1 from LazySets
into format  dot(x,A*x) + 2*dot(b,x) ≤ α
from shape matrix Q and center of ellipsoid c
"""
EllipsoidBBSI(ell::Ellipsoid) = EllipsoidBBSI(ell.center,ell.shape_matrix)


"""
Transform Ellipsoid in format  dot(x,A*x) + 2*dot(b,x) ≤ α to format dot(x-c,Q⁻¹*(x-c)) ≤ 1 from LazySets
from shape matrix Q and center of ellipsoid c
"""
function Ellipsoid(ell::EllipsoidBBSI) 
    c = - (ell.A \ ell.b)
    β = ell.α - dot(c,ell.b)
    Q = Symmetric(inv(Matrix(ell.A)/β))   
   return Ellipsoid(c,Q)
end

"""
This script builds the results and plots presented in XXXX
"""


"""
TestEllipsoids()

"""
function Test_cCRM_TwoEllipsoids(;n :: Int = 100, samples :: Int = 1,
                                  ε :: Real = 1e-6, itmax :: Int = 1000, num_sets :: Int = 2,
                                  restarts = 1, print_file :: Bool = false, 
                                  methods :: Vector{Symbol} = [:centralizedCRM, :MAP], 
                                  useapprox :: Bool = false, 
                                  bench_time :: Bool = false)
    # X = R^n
    # Defines DataFrame for Results
    dfResults, dfFilenames = createDataframes(methods)
    # Fix Random
    Random.seed!(1)
    p = 2*inv(n)
    Random.seed!(1)
    for j in 1:samples
        Ellipsoids = createEllipsoids(n,p,num_sets)
        for i = 1:restarts
            η = -20.
            # x₀ =-20* StartingPoint(n)
            x₀ = fill(η, n)
            prob_name  = savename((Prob=j,Rest=i,n=n,nsets=num_sets))
            timenow= Dates.now()
            dfrow = []
            dfrowFilename = []
            push!(dfrow,prob_name)
            push!(dfrowFilename,prob_name)
            for mtd in methods
                func = eval(mtd) 
                filename = savename("BBIS21",(mtd=mtd,time=timenow),"csv",sort=false)
                print_file ? filedir = datadir("sims",filename) : filedir = ""
                results  = func(x₀,Ellipsoids,itmax=itmax,EPSVAL=ε,gap_distance=true,filedir=filedir)
                elapsed_time = 0.
                if bench_time
                    t = @benchmark $func($x₀,$Ellipsoids,itmax=$itmax,EPSVAL=$ε,gap_distance=true,filedir=$filedir)
                    elapsed_time = (mean(t).time)*1e-9            
                end                
                push!(dfrow,results.iter_total)
                push!(dfrow,elapsed_time)
                push!(dfrowFilename,filedir)
            end    
            push!(dfResults,dfrow)
            push!(dfFilenames,dfrowFilename)
        end
    end
    return dfResults,dfFilenames
end


function Proj_Ellipsoid(x₀::Vector,
                        ell::EllipsoidBBSI)
    
    x₀ ∉ ell ?  nothing : return x₀
    A, b, α  = ell.A, ell.b, ell.α
    n = length(b)
    model = Model()
    @variable(model, x[1:n])
    # set_start_value.(x,x₀)
    @NLexpression(model, Func[i=1:n], x[i] - x₀[i])
    @NLconstraint(model, c, sum(x[i]*sum(A[i,ℓ]*x[ℓ] for  ℓ = 1:n) + 2*b[i]*x[i] for i = 1:n) ≤  α)
    nlp =  MathOptNLSModel(model, Func)
    algencan_specs_file = "algencan.dat"
    writedlm(algencan_specs_file,["ITERATIONS-OUTPUT-DETAIL 00"])
    stats = algencan(nlp,specfnm=algencan_specs_file)
    proj = stats.solution
    return proj
end



function createDataframes(methods::Vector{Symbol})
    dfResults= DataFrame(Problem=String[])
    dfFilenames = copy(dfResults)
    for mtd in methods
        insertcols!(dfResults,join([mtd,"_it"]) => Int[])
        insertcols!(dfResults,join([mtd,"_elapsed"]) => Real[])
        insertcols!(dfFilenames,join([mtd,"filename"]) => String[])
    end
    return dfResults, dfFilenames
end

function createTwoEllipsoids(n::Int, 
                             p::Real; 
                             λ::Real = 1.0)
    Ellipsoids = EllipsoidBBSI[]
    A  = sprandn(n,n,p)
    γ = 1.5
    A = (γ*I + A'*A)
    a = rand(n)
    b = A*a
    adotAa = dot(a,b)
    b .*= -1.
    α = (1+γ)*adotAa
    push!(Ellipsoids,EllipsoidBBSI(A, b, α))
    push!(Ellipsoids, TouchingEllipsoid(Ellipsoids[1], n, λ=λ))
    return Ellipsoids
end

function TouchingEllipsoid(ell::EllipsoidBBSI,
                           n::Int;
                           λ::Real = 1.0)
    c = rand(n)
    while c ∈ ell
        c *= 1.5
    end
    x̂ = Proj_Ellipsoid(c,ell)
    d = λ*(x̂ - c)
    Λ = Diagonal([norm(d); norm(d) .+ 2*rand(n-1)])
    Q, = qr(d)
    A2 = Λ*Q
    M2 = Symmetric(A2'*A2)
    return EllipsoidBBSI(c,M2)
end


function createEllipsoids(n::Int, p::Real, m::Int)
    Ellipsoids = EllipsoidBBSI[]
    for index  in  1:m
        A  = sprandn(n,n,p)
        γ = 1.5
        A = (γ*I + A'*A)
        a = rand(n)
        b = A*a
        adotAa = dot(a,b)
        b .*= -1.
        α = (1+γ)*adotAa
        push!(Ellipsoids,EllipsoidBBSI(A, b, α))
    end
    return Ellipsoids
end


function centralizedCRM(x₀::Vector,
                        Ellipsoids::Vector{EllipsoidBBSI}; 
                        kwargs...)
    ProjectA(x) = Proj_Ellipsoid(x,Ellipsoids[1])
    ProjectB(x) = Proj_Ellipsoid(x,Ellipsoids[2])
    return centralizedCRM(x₀, ProjectA, ProjectB; kwargs...)
end

function MAP(x₀::Vector,
            Ellipsoids::Vector{EllipsoidBBSI}; 
            kwargs...)
    ProjectA(x) = Proj_Ellipsoid(x,Ellipsoids[1])
    ProjectB(x) = Proj_Ellipsoid(x,Ellipsoids[2])
    return MAP(x₀, ProjectA, ProjectB; kwargs...)
end


# # To generate the results.

# size_spaces = [10, 50, 100, 200]
# num_ellipsoids = [5, 10, 20 , 50]
# samples = 10
# # restarts = 1
# ε = 1e-6
# itmax = 50000
# method = [:CRMprod, :MAPprod]
# useapprox = true
# # Too much time. It took 3 days for the results to be finished.
# bench_time = false
# dfResultsEllips, dfEllipFilenames  = createDataframes(method,useapprox)
# for n in size_spaces, m in num_ellipsoids
#     print(m)
#     dfResults, dfFilesname = TestEllipsoids(n=n, num_sets = m, samples = samples, itmax=itmax, 
#                                 ε=ε, bench_time=benhc_time, useapprox=useapprox)
#     append!(dfResultsEllips,dfResults)
#     append!(dfEllipFilenames,dfFilesname)
# end
# ##
# To write data. 
# CSV.write(datadir("sims","AABBIS21_EllipsoidsTable.csv"),dfResultsEllips,
#                     transform = (col, val) -> something(val, missing))
# @show df_describe = describe(dfResultsEllips)
# CSV.write(datadir("sims","AABBIS21_EllipsoidsTableSummary.csv"),
#                     df_describe,transform = (col, val) -> something(val, missing))

# To consult the results.

# dfResultsEllips = CSV.read(datadir("sims","AABBIS21_EllipsoidsTable.csv"), DataFrame)
# @show df_Summary = describe(dfResultsEllips,:mean,:max,:min,:std,:median)[2:end,:]




# perprof = performance_profile(hcat(dfResultsEllips.CRMprodApprox_elapsed, dfResultsEllips.MAPprodApprox_elapsed,
#                                     dfResultsEllips.CRMprod_elapsed,dfResultsEllips.MAPprod_elapsed), 
#                             ["CARM", "MAAP","CRM", "MAP"],
#     title=L"Performance Profile -- Elapsed time comparison -- Gap error -- $\varepsilon = 10^{-6}$",
#     legend = :bottomright, framestyle = :box, linestyles=[:solid, :dash, :dot, :dashdot])
# ylabel!("Percentage of problems solved")
# savefig(perprof,plotsdir("AABBIS21_Ellipsoids_PerProf.pdf"))
# perprof
# perprof2 = performance_profile(hcat(dfResultsEllips.CRMprodApprox_elapsed, dfResultsEllips.MAPprodApprox_elapsed), 
#                             ["CARM", "MAAP"],
#     title=L"Performance Profile -- Elapsed time comparison -- Gap error -- $\varepsilon = 10^{-6}$",
#     legend = :bottomright, framestyle = :box, linestyles=[:solid, :dash,],logscale=false)
# ylabel!("Percentage of problems solved")
# savefig(perprof2,plotsdir("AABBIS21_Ellipsoids_PerProf2.pdf"))
# perprof2
# perprof3 = performance_profile(hcat(dfResultsEllips.CRMprodApprox_elapsed, dfResultsEllips.CRMprod_elapsed), 
#                             ["CARM", "CRM"],
#     title=L"Performance Profile -- Elapsed time comparison -- Gap error -- $\varepsilon = 10^{-6}$",
#     legend = :bottomright, framestyle = :box, linestyles=[:solid, :dash,],logscale=false)
# ylabel!("Percentage of problems solved")
# savefig(perprof3,plotsdir("AABBIS21_Ellipsoids_PerProf3.pdf"))
# perprof3


"""
makeplots_Ellipsoids2D()
Makes plots for Ellipsoids in 2D
"""
function makeplots_Ellipsoids2D(λ::Real, Title::AbstractString)
    ##
    Random.seed!(123)
    n=2
    ℰ = createTwoEllipsoids(n,2/n,λ = λ)
      
    x₀ = [3,2.]
    while x₀ ∈ ℰ[1] || x₀ ∈ ℰ[2]
        x₀ *= 1.2
    end
    resultsMAP = MAP(x₀,ℰ,filedir=datadir("sims/MAP.csv"),itmax = 20000)
    resultscCRM = centralizedCRM(x₀,ℰ,filedir=datadir("sims/centralizedCRM.csv"),itmax = 2000)
    ##
    xMAP = readdlm(datadir("sims/MAP.csv")) 
    xcCRM = readdlm(datadir("sims/centralizedCRM.csv")) 
    plt =  plot(Ellipsoid.(ℰ), framestyle = :none, leg = :bottomleft)
    ##
    plot!([Singleton(v) for v in eachrow(xMAP[2:100,3:4])],label = "MAP -- $(resultsMAP.iter_total) projs.",c = :red,alpha = 0.8, ms = 3)
    plot!([Singleton(v) for v in eachrow(xcCRM[2:end,3:4])],label = "cCRM -- $(resultscCRM.iter_total) projs.", c = :green, alpha = 0.8, ms = 4, m = :diamond)
    MethodPath!(plt,xMAP[1:8,3:4],color = :red)
    MethodPath!(plt,xcCRM[:,3:4],color = :green)
    plot!(Singleton(x₀), c = :black, m = :square, alpha = 1)
    label_points!(plt,[x₀],num_points=1,shift=0.20,label_size=14)
    title!(plt,Title)
    figname = savename("BBIS21_Ellipsoids",(@dict λ),"pdf")
    savefig(plt,plotsdir(figname))
    return plt
end

##
λ = 1.0
Title = L"Ellipsoids intersection without Slater point ($\varepsilon = 10^{-6}$)"
plt1 = makeplots_Ellipsoids2D(λ, Title)
##
λ = 1.1
Title = L"Ellipsoids intersection with Slater point ($\varepsilon = 10^{-6}$)"
plt2 = makeplots_Ellipsoids2D(λ, Title)

