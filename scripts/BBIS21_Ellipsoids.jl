##################################################################
## Basic Functions for Two Ellipsoids tests and plots
##################################################################
using BenchmarkTools
using JuMP, NLPModelsJuMP, NLPModelsAlgencan
import LazySets: Ellipsoid
import Base: in
"""
Structure of an Ellipsoid satisfying dot(x,A*x) + 2*dot(b,x) ≤ α
"""
@kwdef struct EllipsoidBBIS
    A::SparseMatrixCSC
    b::Vector
    α::Float64
end

function in(x₀::Vector, ell::EllipsoidBBIS)
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
function EllipsoidBBIS(c::Vector, Q::AbstractMatrix)
    A = inv(Matrix(Q))
    b = -A*c
    α = 1 + dot(c,b)
    return EllipsoidBBIS(A,b,α)
end


"""
EllipsoidBBIS(ell)

Transform Ellipsoid in format dot(x-c,Q⁻¹*(x-c)) ≤ 1 from LazySets
into format  dot(x,A*x) + 2*dot(b,x) ≤ α
from shape matrix Q and center of ellipsoid c
"""
EllipsoidBBIS(ell::Ellipsoid) = EllipsoidBBIS(ell.center,ell.shape_matrix)


"""
Ellipsoid(ell)
Transform Ellipsoid in format  dot(x,A*x) + 2*dot(b,x) ≤ α to format dot(x-c,Q⁻¹*(x-c)) ≤ 1 from LazySets
from shape matrix Q and center of ellipsoid c
"""
function Ellipsoid(ell::EllipsoidBBIS) 
    c = - (ell.A \ ell.b)
    β = ell.α - dot(c,ell.b)
    Q = Symmetric(inv(Matrix(ell.A)/β))   
   return Ellipsoid(c,Q)
end

"""
Proj_Ellipsoid(x₀, ell)
Projects x₀ onto ell, and EllipsoidBBIS
"""
function Proj_Ellipsoid(x₀::Vector,
                        ell::EllipsoidBBIS)
    x₀ ∉ ell ?  nothing : return x₀
    A, b, α  = ell.A, ell.b, ell.α
    n = length(b)
    model = Model()
    @variable(model, x[1:n])
    @NLexpression(model, Func[i=1:n], x[i] - x₀[i])
    @NLconstraint(model, c, sum(x[i]*sum(A[i,ℓ]*x[ℓ] for  ℓ = 1:n) + 2*b[i]*x[i] for i = 1:n) ≤  α)
    nlp =  MathOptNLSModel(model, Func)
    algencan_specs_file = datadir("exp_pro/algencan.dat")
    writedlm(algencan_specs_file,["ITERATIONS-OUTPUT-DETAIL 00"])
    stats = algencan(nlp,specfnm=algencan_specs_file)
    proj = stats.solution
    return proj
end



function createDataframes(Methods::Vector{Symbol})
    dfResults= DataFrame(Problem=String[])
    dfFilenames = copy(dfResults)
    for mtd in Methods
        insertcols!(dfResults,join([mtd,"_it"]) => Int[])
        insertcols!(dfResults,join([mtd,"_tol"]) => Real[])
        insertcols!(dfResults,join([mtd,"_elapsed"]) => Real[])
        insertcols!(dfFilenames,join([mtd,"filename"]) => String[])
    end
    return dfResults, dfFilenames
end

function createTwoEllipsoids(n::Int, 
                             p::Real; 
                             λ::Real = 1.0)
    Ellipsoids = EllipsoidBBIS[]
    A  = sprandn(n,n,p)
    γ = 1.5
    A = (γ*I + A'*A)
    a = rand(n)
    b = A*a
    adotAa = dot(a,b)
    b .*= -1.
    α = (1+γ)*adotAa
    push!(Ellipsoids,EllipsoidBBIS(A, b, α))
    push!(Ellipsoids, TouchingEllipsoid(Ellipsoids[1], n, λ=λ))
    return Ellipsoids
end

function TouchingEllipsoid(ell::EllipsoidBBIS,
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
    return EllipsoidBBIS(c,M2)
end


function createEllipsoids(n::Int, p::Real, m::Int)
    Ellipsoids = EllipsoidBBIS[]
    for index  in  1:m
        A  = sprandn(n,n,p)
        γ = 1.5
        A = (γ*I + A'*A)
        a = rand(n)
        b = A*a
        adotAa = dot(a,b)
        b .*= -1.
        α = (1+γ)*adotAa
        push!(Ellipsoids,EllipsoidBBIS(A, b, α))
    end
    return Ellipsoids
end

"""
centralizedCRM(x₀, Ellipsoids)
Uses cCRM to find a point into intersection of two EllipsoidBBIS 
"""
function centralizedCRM(x₀::Vector,
                        Ellipsoids::Vector{EllipsoidBBIS}; 
                        kwargs...)
    ProjectA(x) = Proj_Ellipsoid(x,Ellipsoids[1])
    ProjectB(x) = Proj_Ellipsoid(x,Ellipsoids[2])
    return centralizedCRM(x₀, ProjectA, ProjectB; kwargs...)
end


"""
MAP(x₀, Ellipsoids)
Uses MAP to find  a point into intersection of two EllipsoidBBIS 
"""
function MAP(x₀::Vector,
            Ellipsoids::Vector{EllipsoidBBIS}; 
            kwargs...)
    ProjectA(x) = Proj_Ellipsoid(x,Ellipsoids[1])
    ProjectB(x) = Proj_Ellipsoid(x,Ellipsoids[2])
    return MAP(x₀, ProjectA, ProjectB; kwargs...)
end


"""
SPM(x₀, Ellipsoids)
Uses SPM to find  a point into intersection of two EllipsoidBBIS 
"""
function SPM(x₀::Vector,
            Ellipsoids::Vector{EllipsoidBBIS}; 
            kwargs...)
    ProjectA(x) = Proj_Ellipsoid(x,Ellipsoids[1])
    ProjectB(x) = Proj_Ellipsoid(x,Ellipsoids[2])
    return SPM(x₀, ProjectA, ProjectB; kwargs...)
end
##################################################################
## Generating 2D plots
##################################################################
"""
makeplots_Ellipsoids2D()
Makes plots for Ellipsoids in 2D
"""
function makeplots_TwoEllipsoids2D(λ::Real, Title::AbstractString; ε :: Real = 1e-4)
    ##
    Random.seed!(123)
    n=2
    ℰ = createTwoEllipsoids(n,2/n,λ = λ)
    x₀ = [3,2.]
    while x₀ ∈ ℰ[1] || x₀ ∈ ℰ[2]
        x₀ *= 1.2
    end
    resultsMAP = MAP(x₀,ℰ,filedir=datadir("sims/MAP.csv"),itmax = 20_000, EPSVAL=ε)
    resultsSPM = SPM(x₀,ℰ,filedir=datadir("sims/SPM.csv"),itmax = 20_000, EPSVAL=ε)
    resultscCRM = centralizedCRM(x₀,ℰ,filedir=datadir("sims/centralizedCRM.csv"),itmax = 20_000, EPSVAL=ε)
    ##
    xMAP = readdlm(datadir("sims/MAP.csv")) 
    xSPM = readdlm(datadir("sims/SPM.csv")) 
    xcCRM = readdlm(datadir("sims/centralizedCRM.csv")) 
    plt =  plot(Ellipsoid.(ℰ), framestyle = :none, leg = :bottomleft, fillcolor = [:turquoise4  :palegreen4] )
    ##
    plot!([Singleton(v) for v in eachrow(xSPM[2:100,3:4])],label = "SPM -- $(resultsSPM.iter_total) projs.",c = :dodgerblue2,alpha = 0.8, ms = 4, m = :utriangle)
    plot!([Singleton(v) for v in eachrow(xMAP[2:100,3:4])],label = "MAP -- $(resultsMAP.iter_total) projs.",c = :red,alpha = 0.8, ms = 3)
    plot!([Singleton(v) for v in eachrow(xcCRM[2:end,3:4])],label = "cCRM -- $(resultscCRM.iter_total) projs.", c = :green, alpha = 0.8, ms = 4, m = :diamond)
    MethodPath!(plt,xMAP[1:8,3:4], color = :red)
    MethodPath!(plt,xSPM[1:8,3:4], color = :dodgerblue2)
    MethodPath!(plt,xcCRM[:,3:4], color = :green)
    plot!(Singleton(x₀), c = :black, m = :square, alpha = 1)
    label_points!(plt,[x₀],num_points=1,shift=0.20,label_size=14)
    title!(plt,Title)
    figname = savename("BBIS21_Ellipsoids",(@dict λ),"pdf")
    savefig(plt,plotsdir(figname))
    return plt
end
##
    ε = 1e-4
    ε_exponent = floor(Int,log10(ε))
    λ = 1.0
    Title = L"Ellipsoids intersection without Slater point ($\varepsilon = 10^{%$(ε_exponent)}$)"
    plt1 = makeplots_TwoEllipsoids2D(λ, Title, ε = ε)
##
    λ = 1.1
    Title = L"Ellipsoids intersection with Slater point ($\varepsilon = 10^{%$(ε_exponent)}$)"
    plt2 = makeplots_TwoEllipsoids2D(λ, Title, ε = ε)
##################################################################
## Generating Performance Profiles 
###################################################################

"""
TestEllipsoids()

"""
function tests_TwoEllipsoidsRn(;n :: Int = 100, samples :: Int = 1, λ :: Real = 1.0,   
                                  ε :: Real = 1e-6, 
                                  itmax :: Int = 10_000, 
                                  restarts :: Int = 1, 
                                  print_file :: Bool = false, 
                                  Methods :: Vector{Symbol} = [:centralizedCRM, :MAP], 
                                  bench_time :: Bool = false)
    # X = R^n
    # Defines DataFrame for Results
    dfResults, dfFilenames = createDataframes(Methods)
    # Fix Random
    Random.seed!(1)
    # Sparsity of first ellipsoid
    p = 2*inv(n)
    for j in 1:samples
        ℰ = createTwoEllipsoids(n, p, λ = λ)
        for i = 1:restarts
            x₀ = StartingPoint(n)
            while x₀ ∈ ℰ[1] || x₀ ∈ ℰ[2]
                x₀ *= 1.2
            end
            prob_name  = savename((Prob=j,Rest=i,n=n))
            timenow= Dates.now()
            dfrow = []
            dfrowFilename = []
            push!(dfrow,prob_name)
            push!(dfrowFilename,prob_name)
            for mtd in Methods
                func = eval(mtd) 
                filename = savename("BBIS21",(mtd=mtd,time=timenow),"csv",sort=false)
                print_file ? filedir = datadir("sims",filename) : filedir = ""
                results  = func(x₀, ℰ, EPSVAL=ε, gap_distance=true, filedir=filedir, itmax = itmax)
                elapsed_time = 0.
                if bench_time
                    t = @benchmark $func($x₀,$ℰitmax,EPSVAL=$ε,gap_distance=true,filedir=$filedir)
                    elapsed_time = (mean(t).time)*1e-9            
                end                
                push!(dfrow,results.iter_total)
                push!(dfrow,results.final_tol)
                push!(dfrow,elapsed_time)
                push!(dfrowFilename,filedir)
            end    
            push!(dfResults,dfrow)
            push!(dfFilenames,dfrowFilename)
        end
    end
    return dfResults,dfFilenames
end
##
Methods = [:centralizedCRM, :MAP, :SPM]
samples = 30
ε = 1e-4
dfResultsEllips, dfEllipFilenames  = createDataframes(Methods)

## No Slater point in the intersection
dfResults, dfFilesname = tests_TwoEllipsoidsRn(samples = samples, ε = ε, Methods = Methods, λ = 1.0)
append!(dfResultsEllips,dfResults)
append!(dfEllipFilenames,dfFilesname)

## Slater point in the intersection
dfResults, dfFilesname = tests_TwoEllipsoidsRn(samples = samples, ε = ε, Methods = Methods, λ = 1.1)
append!(dfResultsEllips,dfResults)
append!(dfEllipFilenames,dfFilesname)

## To write data. 
CSV.write(datadir("sims","BBIS21_TwoEllipsoidsTable.csv"),dfResultsEllips,
                    transform = (col, val) -> something(val, missing))
@show df_describe = describe(dfResultsEllips)
CSV.write(datadir("sims","BBIS21_TwoEllipsoidsTableSummary.csv"),
                    df_describe,transform = (col, val) -> something(val, missing))

## To consult the results.

dfResultsEllips = CSV.read(datadir("sims","BBIS21_TwoEllipsoidsTable.csv"), DataFrame)
@show df_Summary = describe(dfResultsEllips,:mean,:max,:min,:std,:median)[2:end,:]

## To make Performance Profiles.
perprof = performance_profile(hcat(dfResultsEllips.centralizedCRM_it, dfResultsEllips.MAP_it), 
                            String.(Methods),
                            legend = :bottomright, framestyle = :box, 
                            linestyles=[:solid, :dash, :dot, :dashdot])
ylabel!("Percentage of problems solved")
title!("Performance Profile -- Elapsed time comparison -- Gap error -- $\varepsilon = 10^{-6}$")
savefig(perprof,plotsdir("BBIS21_TwoEllipsoids_Perprof.pdf"))
perprof