__precompile__
using DrWatson
@quickactivate "CRM-CFP"
include(srcdir("CRM-CFP.jl"))
using CSV, DataFrames
using BenchmarkProfiles, LaTeXStrings, Plots
pgfplotsx() 
using BenchmarkTools
using JuMP, NLPModelsJuMP, NLPModelsAlgencan

"""
This script builds the results and plots presented in XXXX
"""


"""
TestEllipsoids()

"""
function TestEllipsoids(;n :: Int = 10, num_sets :: Int = 5, samples :: Int = 1,
        ε :: Real = 1e-6, itmax :: Int = 1000, restarts = 1, print_file :: Bool = false, 
        method :: Vector{Symbol} = [:CRMprod, :MAPprod], useapprox :: Bool = false, bench_time :: Bool = false)
    # X = R^n
    # Defines DataFrame for Results
    dfResults, dfFilenames = createDaframes(method,useapprox)
    # Fix Random
    Random.seed!(1)
    p = 2*inv(n)
    Random.seed!(1)
    for j in 1:samples
        Ellipsoids = createEllipsoids(n,p,num_sets)
        for i = 1:restarts
            η = -2.
            # x₀ = StartingPoint(n)
            x₀ = fill(η, n)
            prob_name  = savename((Prob=j,Rest=i,n=n,nsets=num_sets))
            timenow= Dates.now()
            dfrow = []
            dfrowFilename = []
            push!(dfrow,prob_name)
            push!(dfrowFilename,prob_name)
            for mtd in method
                func = eval(mtd) 
                filename = savename("AABBIS21",(mtd=mtd,time=timenow),"csv",sort=false)
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
                if useapprox
                    Ellipsoids[begin][:useapprox] = true
                    mtd = Symbol("Approx"*String(mtd))
                    filename = savename("AABBIS21",(mtd=mtd,time=timenow),"csv",sort=false)
                    print_file ? filedir = datadir("sims",filename) : filedir = ""
                    results  = func(x₀,Ellipsoids,itmax=itmax,EPSVAL=ε,gap_distance=true,filedir=filedir)
                    elapsed_time = 0.0
                    if bench_time
                        t = @benchmark $func($x₀,$Ellipsoids,itmax=$itmax,EPSVAL=$ε,gap_distance=true,filedir=$filedir)
                        elapsed_time = (mean(t).time)*1e-9            
                    end
                    push!(dfrow,results.iter_total)
                    push!(dfrow,elapsed_time)
                    push!(dfrowFilename,filedir)
                    Ellipsoids[begin][:useapprox] = false

                end
            end    
            push!(dfResults,dfrow)
            push!(dfFilenames,dfrowFilename)
        end
    end
    return dfResults,dfFilenames
end


function ProjectProdSpace(X::Vector,Ellipsoids::Vector{Dict})
    proj = similar(X)
    useapprox = Ellipsoids[begin][:useapprox]
    for index in eachindex(proj)
        func = useapprox ? eval(:ApproxProj_Ellipsoid) : eval(:Proj_Ellipsoid)
        proj[index] = func(X[index],Ellipsoids[index])
    end
    return proj
end


# function ProjectProdSpace(X::Vector,Ellipsoids::Vector{Dict})
#     proj = similar(X)
#     for index in eachindex(proj)
#         proj[index] = ApproxProj_Ellipsoid(X[index],Ellipsoids[index])
#     end
#     return proj
# end


function  ApproxProj_Ellipsoid(x::Vector,Ellipsoid::Dict; λ::Real = 1.0)
  @unpack A, b, α  = Ellipsoid
  Ax = A*x
  gx = dot(x,Ax) + 2*dot(b,x) - α
  if gx ≤ 0
    return x
  else
    ∂gx = 2*(Ax + b)
    return λ*( x .- (gx/dot(∂gx,∂gx))*∂gx ).+ (1-λ)*x
  end
end

function Proj_Ellipsoid(x₀::Vector,Ellipsoid::Dict)
  @unpack A, b, α  = Ellipsoid
  if dot(x₀,A*x₀) + 2*dot(b,x₀) ≤ α
    return x₀
  else
    n = length(b)
    model = Model()
    @variable(model, x[1:n])
    # set_start_value(x,x₀)
    @NLexpression(model, Func[i=1:n], x[i] - x₀[i])
    @NLconstraint(model, c, sum(x[i]*sum(A[i,ℓ]*x[ℓ] for  ℓ=1:n) + 2*b[i]*x[i] for i = 1:n) ≤  α)
    nlp =  MathOptNLSModel(model, Func)
    algencan_specs_file = "algencan.dat"
    writedlm(algencan_specs_file,["ITERATIONS-OUTPUT-DETAIL 00"])
    stats = algencan(nlp,specfnm=algencan_specs_file)
    proj = stats.solution
    return proj
  end
end


function createDaframes(method::Vector{Symbol},useapprox::Bool)
    dfResults= DataFrame(Problem=String[])
    dfFilenames = copy(dfResults)
    for mtd in method
        insertcols!(dfResults,join([mtd,"_it"]) => Int[])
        insertcols!(dfResults,join([mtd,"_elapsed"]) => Real[])
        insertcols!(dfFilenames,join([mtd,"filename"]) => String[])
        if useapprox
            insertcols!(dfResults,join([mtd,"Approx_it"]) => Int[])
            insertcols!(dfResults,join([mtd,"Approx_elapsed"]) => Real[])
            insertcols!(dfFilenames,join([mtd,"Approxfilename"]) => String[])
        end
    end
    return dfResults, dfFilenames
end

function createEllipsoids(n::Int, p::Real, m::Int)
    Ellipsoids = Dict[]
    for index  in  1:m
        A  = sprandn(n,n,p)
        γ = 1.5
        A = (γ*I + A'*A)
        a = rand(n)
        b = A*a
        adotAa = dot(a,b)
        b .*= -1.
        α = (1+γ)*adotAa
        push!(Ellipsoids,Dict([:A, :b, :α, :useapprox] .=> [A, b, α, false]))
    end
    return Ellipsoids
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
# dfResultsEllips, dfEllipFilenames  = createDaframes(method,useapprox)
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

dfResultsEllips = CSV.read(datadir("sims","AABBIS21_EllipsoidsTable.csv"), DataFrame)
@show df_Summary = describe(dfResultsEllips,:mean,:max,:min,:std,:median)[2:end,:]




perprof = performance_profile(hcat(dfResultsEllips.CRMprodApprox_elapsed, dfResultsEllips.MAPprodApprox_elapsed,
                                    dfResultsEllips.CRMprod_elapsed,dfResultsEllips.MAPprod_elapsed), 
                            ["CARM", "MAAP","CRM", "MAP"],
    title=L"Performance Profile -- Elapsed time comparison -- Gap error -- $\varepsilon = 10^{-6}$",
    legend = :bottomright, framestyle = :box, linestyles=[:solid, :dash, :dot, :dashdot])
ylabel!("Percentage of problems solved")
savefig(perprof,plotsdir("AABBIS21_Ellipsoids_PerProf.pdf"))
perprof
perprof2 = performance_profile(hcat(dfResultsEllips.CRMprodApprox_elapsed, dfResultsEllips.MAPprodApprox_elapsed), 
                            ["CARM", "MAAP"],
    title=L"Performance Profile -- Elapsed time comparison -- Gap error -- $\varepsilon = 10^{-6}$",
    legend = :bottomright, framestyle = :box, linestyles=[:solid, :dash,],logscale=false)
ylabel!("Percentage of problems solved")
savefig(perprof2,plotsdir("AABBIS21_Ellipsoids_PerProf2.pdf"))
perprof2
perprof3 = performance_profile(hcat(dfResultsEllips.CRMprodApprox_elapsed, dfResultsEllips.CRMprod_elapsed), 
                            ["CARM", "CRM"],
    title=L"Performance Profile -- Elapsed time comparison -- Gap error -- $\varepsilon = 10^{-6}$",
    legend = :bottomright, framestyle = :box, linestyles=[:solid, :dash,],logscale=false)
ylabel!("Percentage of problems solved")
savefig(perprof3,plotsdir("AABBIS21_Ellipsoids_PerProf3.pdf"))
perprof3







