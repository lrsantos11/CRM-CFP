using BenchmarkTools
using JuMP, NLPModelsJuMP, NLPModelsAlgencan


"""
Structure of an Ellipsoid satisfying dot(x,A*x) + 2*dot(b,x) ≤ α
"""
@kwdef struct EllipsoidI
    A::SparseMatrixCSC
    b::Vector
    α::Float64
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
                        ell::EllipsoidI)
  A, b, α  = ell.A, ell.b, ell.α
  if dot(x₀,A*x₀) + 2*dot(b,x₀) ≤ α
    return x₀
  else
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

function createEllipsoids(n::Int, p::Real, m::Int)
    Ellipsoids = EllipsoidI[]
    for index  in  1:m
        A  = sprandn(n,n,p)
        γ = 1.5
        A = (γ*I + A'*A)
        a = rand(n)
        b = A*a
        adotAa = dot(a,b)
        b .*= -1.
        α = (1+γ)*adotAa
        push!(Ellipsoids,EllipsoidI(A, b, α))
    end
    return Ellipsoids
end


function centralizedCRM(x₀::Vector,
                        Ellipsoids::Vector{EllipsoidI}; 
                        kwargs...)
    ProjectA(x) = Proj_Ellipsoid(x,Ellipsoids[1])
    ProjectB(x) = Proj_Ellipsoid(x,Ellipsoids[2])
    return centralizedCRM(x₀, ProjectA, ProjectB; kwargs...)
end

function MAP(x₀::Vector,
            Ellipsoids::Vector{EllipsoidI}; 
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
function makeplots_Ellipsoids2D()
    ##
    ##### Repeating Ellipsoid Examples of the paper  Approximate Douglas-Rachford algorithm for two-sets convex feasibility problems, by R. Díaz Millán, O.P. Ferreira, J. Ugon (arXiv:2105.13005)

    ##
    # Rotation Matrix
    Rot(θ) = [cos(θ) sin(θ);-sin(θ) cos(θ)]

    # First ellipse. We will find the intersection between this set and 4 other sets.

    A = diagm([2.0,1/5.0])*Rot(-π/4)
    M = inv(Symmetric(A'*A))
    z0=[0.0, 0.]
    ell = LazySets.Ellipsoid(z0,inv(M)) 
    ell_BBIS = EllipsoidI(M,z0, 1)


    ##
    # One ellipse that intersects with the first ellipse, with Slater point.
    B = diagm([1/2.0,5.0/2])*Rot(π/3)
    M1 = Symmetric(B'*B)
    z1 = [2.3,1.0/2]
    ell1 = LazySets.Ellipsoid(z1,inv(M1))
    ell1_BBIS = EllipsoidI((M1),-M1*z1, 1 - dot(z1,M1*z1))
    Ellipsoids_Ex1 = [ell_BBIS; ell1_BBIS]
    x₀ = [-1,1.5]
    resultsMAP = MAP(x₀,Ellipsoids_Ex1,filedir="MAP.csv",print_intermediate=true, itmax = 4000, gap_distance = true, EPSVAL = 1e-6)
    resultscCRM = centralizedCRM(x₀,Ellipsoids_Ex1,filedir="cCRM.csv", itmax = 400, gap_distance = true, EPSVAL = 1e-6)
    pointsMAP  = readdlm("MAP.csv")
    pointscCRM  = readdlm("cCRM.csv")
    ##
    plt1 = plot(ell,framestyle = :none, aspect_ratio = :equal , leg = :bottomleft)
    plot!(plt1,ell1)
    plot!(plt1,[Singleton(v) for v in eachrow(pointsMAP[:,3:4])], label="MAP -- $(resultsMAP.iter_total) proj")
    plot!(plt1,[Singleton(v) for v in eachrow(pointscCRM[:,3:4])],c=:red,label="cCRM -- $(resultscCRM.iter_total) proj")
    MethodPath!(plt1,pointsMAP[:,3:4],color = :green)
    MethodPath!(plt1,pointscCRM[:,3:4],color = :red)
    title!(plt1,L"Ellipsoids with  Slater point in intersection ($\varepsilon = 10^{-6}$)")
    label_points!(plt1,[x₀],num_points=1,label_size=12,shift=0.15)
    plot!(plt1,Singleton(x₀), m=:square,c=:black,alpha=1)



    ##
    # One ellipse that intersects with the first ellipse, without Slater point.
    d = [1.273,-1.1968]         # Pick a direction
    x = d/sqrt(dot((M)*d,d)) # Find a point on the boundary
    dd = inv(M1)*(M)*x 
    λ = dd'*M1*dd
    z2 = dd/λ + x
    M2 = M1*λ
    ell2 = LazySets.Ellipsoid(z2,inv(M2))
    ell2_BBIS = EllipsoidI((M2),-M2*z2, 1 - dot(z2,M2*z2))
    plt2 = plot(ell,framestyle = :none, aspect_ratio = :equal , leg = :bottomleft)
    plot!(plt2,ell2)
    title!(plt2,L"Ellipsoids without  Slater point in intersection ($\varepsilon = 10^{-6}$)")

    ##
    Ellipsoids_Ex2 = [ell_BBIS; ell2_BBIS]
    x₀ = [-1,1.5]
    resultsMAP = MAP(x₀,Ellipsoids_Ex2,filedir="MAP.csv",print_intermediate=true, itmax = 100_000, gap_distance = true, EPSVAL = 1e-6)
    resultscCRM = centralizedCRM(x₀,Ellipsoids_Ex2,filedir="cCRM.csv", itmax = 400, gap_distance = true, EPSVAL = 1e-6)
    ##
    pointsMAP  = readdlm("MAP.csv")
    pointscCRM  = readdlm("cCRM.csv")
    plot!(plt2,[Singleton(v) for v in eachrow(pointsMAP[1:30,3:4])], label="MAP -- $(resultsMAP.iter_total) proj")
    plot!(plt2,[Singleton(v) for v in eachrow(pointscCRM[:,3:4])],c=:red,label="cCRM -- $(resultscCRM.iter_total) proj")
    MethodPath!(plt2,pointsMAP[1:20,3:4],color = :green)
    MethodPath!(plt2,pointscCRM[:,3:4],color = :red)
    label_points!(plt1,[x₀],num_points=1,label_size=12,shift=0.15)
    plot!(plt1,Singleton(x₀), m=:square,c=:black,alpha=1)

    ##
    savefig(plt1,plotsdir("BBIS21_Ellipsoids_withSlaterPoint.pdf"))
    savefig(plt2,plotsdir("BBIS21_Ellipsoids_withhtouSlaterPoint.pdf"))
end

