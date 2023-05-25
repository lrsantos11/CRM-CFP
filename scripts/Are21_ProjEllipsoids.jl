__precompile__
using DrWatson
@quickactivate "CRM-CFP"
include(srcdir("CRM-CFP.jl"))
using BenchmarkTools, BenchmarkProfiles
using CSV, DataFrames
using LaTeXStrings, Plots
pgfplotsx()
import LazySets: Ellipsoid
import Base: in
##################################################################
## Basic Functions for Ellipsoids tests and plots
##################################################################





function ProjectEllipsoids_ProdSpace(X::Vector,
    Ellipsoids::Vector{EllipsoidCRM})
    proj = similar(X)
    for index in eachindex(proj)
        proj[index] = Proj_Ellipsoid(X[index], Ellipsoids[index])
    end
    return proj
end


function ApproxProjectEllipsoids_ProdSpace(X::Vector,
    Ellipsoids::Vector{EllipsoidCRM})
    proj = similar(X)
    for index in eachindex(proj)
        proj[index] = ApproxProj_Ellipsoid(X[index], Ellipsoids[index])
    end
    return proj
end




function ApproxProj_Ellipsoid(x::Vector,
    Ellipsoid::EllipsoidCRM;
    λ::Real = 1.0)
    A, b, α = Ellipsoid.A, Ellipsoid.b, Ellipsoid.α
    Ax = A * x
    gx = dot(x, Ax) + 2 * dot(b, x) - α
    if gx ≤ 0
        return x
    else
        ∂gx = 2 * (Ax + b)
        return λ * (x .- (gx / dot(∂gx, ∂gx)) * ∂gx) .+ (1 - λ) * x
    end
end


function createEllipsoids(n::Int, p::Real, num_sets::Int)
    Ellipsoids = EllipsoidCRM[]
    for index = 1:num_sets
        A = sprandn(n, n, p)
        γ = 1.5
        A = (γ * I + A' * A)
        a = rand(n)
        b = A * a
        adotAa = dot(a, b)
        b .*= -1.0
        α = (1 + γ) * adotAa
        push!(Ellipsoids, EllipsoidCRM(A, b, α))
    end
    return Ellipsoids
end
"""
make_blocks()

"""

function make_blocks(num_sets::Int)
    max_size = 5
    min_size = 3
    block_indexes = Vector[]
    flag = true
    indexes = collect(1:num_sets)
    sum_blocks = 0
    while sum_blocks != num_sets && flag
        blocksize = rand(min_size:max_size)
        sum_blocks += blocksize
        if sum_blocks > num_sets
            blocksize -= (sum_blocks - num_sets)
            flag = false
        end
        push!(block_indexes, indexes[1:blocksize])
        deleteat!(indexes, 1:blocksize)
    end
    return block_indexes
end


function generate_fneProj(block_indexes, num_sets, Ellipsoids)
    if isempty(block_indexes)
        block_indexes = make_blocks(num_sets)
    end
    num_blocks = length(block_indexes)
    Projections = Function[]
    for block in block_indexes
        push!(Projections, x -> mean([Proj_Ellipsoid(x, ell) for ell in Ellipsoids[block]]))
    end
    return Projections
end


##################################################################
## Methods specialized for this set of experiments.
##################################################################


"""
CRMprod(x₀, Ellipsoids)
Uses CRMprod to find a point into intersection of a finite number of  EllipsoidCRM 
"""
function CRMprod(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM};
    block_indexes::Vector{Vector} = Vector[],
    kwargs...)
    Projections = Function[]
    for ell in Ellipsoids
        push!(Projections, x -> Proj_Ellipsoid(x, ell))
    end
    return CRMprod(x₀, Projections; kwargs...)
end


"""
fneCRMprod(x₀, Ellipsoids)
Uses fneCRMprod to find a point into intersection of a finite number of  EllipsoidCRM 
"""
function fneCRMprod(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM};
    block_indexes::Vector{Vector} = Vector[],
    kwargs...)
    num_sets = length(Ellipsoids)
    if isempty(block_indexes)
        block_indexes = make_blocks(num_sets)
    end
    Projections = generate_fneProj(block_indexes, num_sets, Ellipsoids)
    return CRMprod(x₀, Projections; kwargs...)
end



"""
fneMAPprod(x₀, Ellipsoids)
Uses fneMAPprod to find a point into intersection of a finite number of  EllipsoidCRM 
"""
function fneCimmino(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM};
    block_indexes::Vector{Vector} = Vector[],
    kwargs...)
    num_sets = length(Ellipsoids)
    if isempty(block_indexes)
        block_indexes = make_blocks(num_sets)
    end
    Projections = generate_fneProj(block_indexes, num_sets, Ellipsoids)
    return MAPprod(x₀, Projections; kwargs...)
end

"""
CARMprod(x₀, Ellipsoids)
Uses CARMprod to find a point into intersection of a finite number of  EllipsoidCRM 
"""
function CARMprod(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM};
    block_indexes::Vector{Vector} = Vector[],
    kwargs...)
    Projections = Function[]
    for ell in Ellipsoids
        push!(Projections, x -> ApproxProj_Ellipsoid(x, ell))
    end
    return CRMprod(x₀, Projections; kwargs...)
end




"""
MAP(x₀, Ellipsoids)
Uses MAP to find  a point into intersection of a finite number of  EllipsoidCRM 
"""
function MAPprod(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM};
    block_indexes::Vector{Vector} = Vector[],
    kwargs...)
    Projections = Function[]
    for ell in Ellipsoids
        push!(Projections, x -> Proj_Ellipsoid(x, ell))
    end
    return MAPprod(x₀, Projections; kwargs...)
end



"""
MAAP(x₀, Ellipsoids)
Uses MAAP to find  a point into intersection of a finite number of  EllipsoidCRM 
"""
function MAAPprod(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM};
    block_indexes::Vector{Vector} = Vector[],
    kwargs...)
    Projections = Function[]
    for ell in Ellipsoids
        push!(Projections, x -> ApproxProj_Ellipsoid(x, ell))
    end
    return MAPprod(x₀, Projections; kwargs...)
end

##################################################################
## Functions to run the benchmark
##################################################################


"""
createDataFrames()

"""
function createDataFrames(method::Vector{Symbol}, bench_time::Bool)
    dfResults = DataFrame(Problem = String[])
    dfFilenames = copy(dfResults)
    insertcols!(dfResults, :OpBlocks => Vector[])
    insertcols!(dfResults, :Num_Operators => Int[])
    for mtd in method
        insertcols!(dfResults, join([mtd, "_it"]) => Int[])
        bench_time ? insertcols!(dfResults, join([mtd, "_elapsed"]) => Real[]) : nothing
        insertcols!(dfFilenames, join([mtd, "filename"]) => String[])
    end
    return dfResults, dfFilenames
end

"""
Benchmark_Ellipsoids()

"""
function Benchmark_Ellipsoids(;
    n::Int = 100,
    num_sets::Int = 10,
    samples::Int = 1,
    ε::Real = 1e-6,
    itmax::Int = 1000,
    restarts::Int = 1,
    print_file::Bool = false,
    method::Vector{Symbol} = [:CRMprod, :fneCRMprod, :fneCimmino],
    bench_time::Bool = false)
    # X = R^n
    # Defines DataFrame for Results
    dfResults, dfFilenames = createDataFrames(method, bench_time)
    # Fix Random
    Random.seed!(1)
    p = 2 * inv(n)
    for j in 1:samples
        Ellipsoids = createEllipsoids(n, p, num_sets)
        for i = 1:restarts
            η = -2.0
            # x₀ = StartingPoint(n)
            x₀ = fill(η, n)
            prob_name = savename((Prob = j, Rest = i, n = n, nsets = num_sets))
            timenow = Dates.now()
            dfrow = []
            dfrowFilename = []
            push!(dfrow, prob_name)
            push!(dfrowFilename, prob_name)
            block_indexes = make_blocks(num_sets)
            push!(dfrow, block_indexes)
            num_op = length(block_indexes)
            push!(dfrow, num_op)
            for mtd in method
                func = eval(mtd)
                filename = savename("Are21", (mtd = mtd, time = timenow), "csv", sort = false)
                print_file ? filedir = datadir("sims", filename) : filedir = ""
                results = func(x₀, Ellipsoids, itmax = itmax, EPSVAL = ε, gap_distance = true, filedir = filedir, block_indexes = block_indexes)
                push!(dfrow, results.iter_total)
                if bench_time
                    t = @benchmark $func($x₀, $Ellipsoids, itmax = $itmax, EPSVAL = $ε, gap_distance = true, filedir = $filedir, block_indexes = $block_indexes)
                    elapsed_time = (mean(t).time) * 1e-9
                    push!(dfrow, elapsed_time)

                end
                push!(dfrowFilename, filedir)
            end
            push!(dfResults, dfrow)
            push!(dfFilenames, dfrowFilename)
        end
    end
    return dfResults, dfFilenames
end

"""
Benchmark_Ellipsoids()

"""
function Benchmark_Ellipsoids_fne(;
    n::Int = 100,
    num_sets::Int = 10,
    samples::Int = 1,
    ε::Real = 1e-6,
    itmax::Int = 10000,
    restarts::Int = 5,
    print_file::Bool = false,
    method::Vector{Symbol} = [:fneCRMprod, :fneCimmino],
    bench_time::Bool = false)
    # X = R^n
    # Defines DataFrame for Results
    dfResults, dfFilenames = createDataFrames(method, bench_time)
    # Fix Random
    Random.seed!(1)
    p = 2 * inv(n)
    for j in 1:samples
        block_indexes = make_blocks(num_sets)
        num_op = length(block_indexes)
        for i = 1:restarts
            Ellipsoids = createEllipsoids(n, p, num_sets)
            η = -2.0
            # x₀ = StartingPoint(n)
            x₀ = fill(η, n)
            prob_name = savename((Block = j, Prob = i, n = n, nsets = num_sets))
            timenow = Dates.now()
            dfrow = []
            dfrowFilename = []
            push!(dfrow, prob_name)
            push!(dfrowFilename, prob_name)
            push!(dfrow, block_indexes)
            push!(dfrow, num_op)
            for mtd in method
                func = eval(mtd)
                filename = savename("Are21", (mtd = mtd, time = timenow), "csv", sort = false)
                print_file ? filedir = datadir("sims", filename) : filedir = ""
                results = func(x₀, Ellipsoids, itmax = itmax, EPSVAL = ε, gap_distance = true, filedir = filedir, block_indexes = block_indexes)
                push!(dfrow, results.iter_total)
                if bench_time
                    t = @benchmark $func($x₀, $Ellipsoids, itmax = $itmax, EPSVAL = $ε, gap_distance = true, filedir = $filedir, block_indexes = $block_indexes)
                    elapsed_time = (mean(t).time) * 1e-9
                    push!(dfrow, elapsed_time)

                end
                push!(dfrowFilename, filedir)
            end
            push!(dfResults, dfrow)
            push!(dfFilenames, dfrowFilename)
        end
    end
    return dfResults, dfFilenames
end

##################################################################
## Benchmark #1 - CRMProd vs MAPprod vs fneCRMprod
##################################################################


# # To generate the results.

size_spaces = [10, 30, 50, 100, 200]
num_ellipsoids = [10, 30, 50, 100, 200]
samples = 10
# restarts = 1
ε = 1e-6
itmax = 50000
methods = [:CRMprod, :fneCRMprod, :fneCimmino, :MAPprod, :CARMprod, :MAAPprod]
# # Too much time. It took 1 day for the results to be finished.
bench_time = true
dfResultsEllips, dfEllipFilenames = createDataFrames(methods, bench_time)
#for n in size_spaces, num_sets in num_ellipsoids
 #   println("Running benchmark on $num_sets ellipsoids in R^$(n)")
  #  dfResults, dfFilesname = Benchmark_Ellipsoids(n = n, num_sets = num_sets, samples = samples, itmax = itmax, ε = ε, bench_time = bench_time, method = methods)
  #  append!(dfResultsEllips, dfResults)
  #  append!(dfEllipFilenames, dfFilesname)
#end

##
# Write data. 
timenow = Dates.now()
CSV.write(datadir("sims", savename("Are21_EllipsoidsTable", (time = timenow,), "csv", sort = false)), dfResultsEllips, transform = (col, val) -> something(val, missing))
@show df_describe = describe(dfResultsEllips, :min, :median, :mean, :max, :std)
CSV.write(datadir("sims", savename("Are21_EllipsoidsTableSummary", (time = timenow,), "csv", sort = false)),
    df_describe, transform = (col, val) -> something(val, missing))




##################################################################
## Benchmark #2 - fneCRMprod vs fneCimmino - number of 
##################################################################


# # To generate the results.

size_spaces = [10, 30, 50, 100, 200]
num_ellipsoids = [30, 50, 100, 200]
samples_blocks = 3
restarts = 10
ε = 1e-6
itmax = 50000
methods = [:fneCRMprod, :fneCimmino]
# # # Too much time. It took 3 days for the results to be finished.
bench_time = true
dfResultsEllips, dfEllipFilenames = createDataFrames(methods, bench_time)
for n in size_spaces, num_sets in num_ellipsoids
    println("Running benchmark on $num_sets ellipsoids in R^$(n)")
    dfResults, dfFilesname = Benchmark_Ellipsoids_fne(n = n, num_sets = num_sets, samples = samples_blocks, itmax = itmax, ε = ε, bench_time = bench_time, method = methods, restarts = restarts)
    append!(dfResultsEllips, dfResults)
    append!(dfEllipFilenames, dfFilesname)
end


##
# Write data. 
timenow = Dates.now()
CSV.write(datadir("sims", savename("Are21_Ellipsoids_fneTest", (time = timenow,), "csv", sort = false)), dfResultsEllips, transform = (col, val) -> something(val, missing))
@show df_describe = describe(dfResultsEllips, :min, :median, :mean, :max, :std)
CSV.write(datadir("sims", savename("Are21_Ellipsoids_fneTestTableSummary", (time = timenow,), "csv", sort = false)),
    df_describe, transform = (col, val) -> something(val, missing))




# To consult the results.

# dfResultsEllips = CSV.read(datadir("sims","Are21_EllipsoidsTable.csv"), DataFrame)
# @show df_Summary = describe(dfResultsEllips,:mean,:max,:min,:std,:median)[2:end,:]



## Generate Peformance Profiles

# perprof = performance_profile(hcat(dfResultsEllips.CRMprodApprox_elapsed, dfResultsEllips.MAPprodApprox_elapsed,
#                                     dfResultsEllips.CRMprod_elapsed,dfResultsEllips.MAPprod_elapsed), 
#                             ["CARM", "MAAP","CRM", "MAP"],
#     title=L"Performance Profile -- Elapsed time comparison -- Gap error -- $\varepsilon = 10^{-6}$",
#     legend = :bottomright, framestyle = :box, linestyles=[:solid, :dash, :dot, :dashdot])
# ylabel!("Percentage of problems solved")
# savefig(perprof,plotsdir("Are21_Ellipsoids_PerProf.pdf"))
# perprof
# perprof2 = performance_profile(hcat(dfResultsEllips.CRMprodApprox_elapsed, dfResultsEllips.MAPprodApprox_elapsed), 
#                             ["CARM", "MAAP"],
#     title=L"Performance Profile -- Elapsed time comparison -- Gap error -- $\varepsilon = 10^{-6}$",
#     legend = :bottomright, framestyle = :box, linestyles=[:solid, :dash,],logscale=false)
# ylabel!("Percentage of problems solved")
# savefig(perprof2,plotsdir("Are21_Ellipsoids_PerProf2.pdf"))
# perprof2
# perprof3 = performance_profile(hcat(dfResultsEllips.CRMprodApprox_elapsed, dfResultsEllips.CRMprod_elapsed), 
#                             ["CARM", "CRM"],
#     title=L"Performance Profile -- Elapsed time comparison -- Gap error -- $\varepsilon = 10^{-6}$",
#     legend = :bottomright, framestyle = :box, linestyles=[:solid, :dash,],logscale=false)
# ylabel!("Percentage of problems solved")
# savefig(perprof3,plotsdir("Are21_Ellipsoids_PerProf3.pdf"))
# perprof3



##


