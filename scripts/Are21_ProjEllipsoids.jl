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

"""
Structure of an Ellipsoid satisfying dot(x,A*x) + 2*dot(b,x) ≤ α
"""
@kwdef struct EllipsoidBBIS
    A::AbstractMatrix
    b::AbstractVector
    α::Number
end

function in(x₀::Vector, ell::EllipsoidBBIS)
    A, b, α = ell.A, ell.b, ell.α
    if dot(x₀, A * x₀) + 2 * dot(b, x₀) ≤ α
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
    b = -A * c
    α = 1 + dot(c, b)
    return EllipsoidBBIS(A, b, α)
end


"""
EllipsoidBBIS(ell)

Transform Ellipsoid in format dot(x-c,Q⁻¹*(x-c)) ≤ 1 from LazySets
into format  dot(x,A*x) + 2*dot(b,x) ≤ α
from shape matrix Q and center of ellipsoid c
"""
EllipsoidBBIS(ell::Ellipsoid) = EllipsoidBBIS(ell.center, ell.shape_matrix)


"""
Ellipsoid(ell)
Transform Ellipsoid in format  dot(x,A*x) + 2*dot(b,x) ≤ α to format dot(x-c,Q⁻¹*(x-c)) ≤ 1 from LazySets
from shape matrix Q and center of ellipsoid c
"""
function Ellipsoid(ell::EllipsoidBBIS)
    c = -(ell.A \ ell.b)
    β = ell.α - dot(c, ell.b)
    Q = Symmetric(inv(Matrix(ell.A) / β))
    return Ellipsoid(c, Q)
end



"""
Proj_Ellipsoid(x₀, ell)
Projects x₀ onto ell, and EllipsoidBBIS using an ADMM algorithm as reported by Jia, Cai and Han [Jia2007]

[Jia2007] Z. Jia, X. Cai, e D. Han, “Comparison of several fast algorithms for projection onto an ellipsoid”, Journal of Computational and Applied Mathematics, vol. 319, p. 320–337, ago. 2017, doi: 10.1016/j.cam.2017.01.008.
"""
function Proj_Ellipsoid(x₀::Vector,
    ell::EllipsoidBBIS;
    itmax::Int = 10_000,
    ε::Real = 1e-8)
    x₀ ∉ ell ? nothing : return x₀
    A, b, α = ell.A, ell.b, ell.α
    ϑₖ = 10 / norm(A)
    n = length(x₀)
    B = sqrt(Matrix(A))
    issymmetric(B) ? BT = B : BT = B'
    b̄ = B \ (-b)
    αplusb̄2 = α + norm(b̄)^2
    r = sqrt(αplusb̄2)
    yₖ = ones(n)
    λₖ = ones(n)
    xₖ = x₀
    it = 0
    tolADMM = 1.0
    function ProjY(y)
        normy = norm(y)
        if αplusb̄2 - normy^2 ≥ 0.0
            return y
        else
            return (r / normy) * y
        end
    end
    Ā = (I + ϑₖ * A)
    normRxₖ = 0.0
    while tolADMM ≥ ε^2 && it ≤ itmax
        uₖ = x₀ + BT * (λₖ + ϑₖ * (yₖ + b̄))
        xₖ = Ā \ uₖ
        wₖ = B * xₖ - λₖ ./ ϑₖ - b̄
        normwₖ = norm(wₖ)
        normwₖ ≤ r ? yₖ = wₖ : yₖ = (r / normwₖ) * wₖ
        Rxₖ = xₖ - x₀ - BT * λₖ
        Ryₖ = yₖ - ProjY(yₖ - λₖ)
        Rλₖ = B * xₖ - yₖ - b̄
        λₖ -= ϑₖ * Rλₖ
        normRxₖ = norm(Rxₖ)
        normRyₖ = norm(Ryₖ)
        normRλₖ = norm(Rλₖ)
        tolADMM = sum([normRxₖ^2, normRyₖ^2, normRλₖ^2])
        it += 1
        # if normRxₖ < normRλₖ*(0.1/n)
        #     ϑₖ *= 2
        # elseif normRxₖ > normRλₖ*(0.9/n)
        #     ϑₖ *= 0.5
        # end

    end

    # @show it, normRxₖ
    return real.(xₖ)
end





function ProjectEllipsoids_ProdSpace(X::Vector,
                                     Ellipsoids::Vector{EllipsoidBBIS})
    proj = similar(X)
    for index in eachindex(proj)
        proj[index] = Proj_Ellipsoid(X[index], Ellipsoids[index])
    end
    return proj
end


function ApproxProjectEllipsoids_ProdSpace(X::Vector,
                                           Ellipsoids::Vector{EllipsoidBBIS})
    proj = similar(X)
    for index in eachindex(proj)
        proj[index] = ApproxProj_Ellipsoid(X[index], Ellipsoids[index])
    end
    return proj
end




function ApproxProj_Ellipsoid(x::Vector,
    Ellipsoid::EllipsoidBBIS;
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


function createDaframes(method::Vector{Symbol}, useapprox::Bool)
    dfResults = DataFrame(Problem = String[])
    dfFilenames = copy(dfResults)
    for mtd in method
        insertcols!(dfResults, join([mtd, "_it"]) => Int[])
        insertcols!(dfResults, join([mtd, "_elapsed"]) => Real[])
        insertcols!(dfFilenames, join([mtd, "filename"]) => String[])
        if useapprox
            insertcols!(dfResults, join([mtd, "Approx_it"]) => Int[])
            insertcols!(dfResults, join([mtd, "Approx_elapsed"]) => Real[])
            insertcols!(dfFilenames, join([mtd, "Approxfilename"]) => String[])
        end
    end
    return dfResults, dfFilenames
end

function createEllipsoids(n::Int, p::Real, num_sets::Int)
    Ellipsoids = EllipsoidBBIS[]
    for index = 1:num_sets
        A = sprandn(n, n, p)
        γ = 1.5
        A = (γ * I + A' * A)
        a = rand(n)
        b = A * a
        adotAa = dot(a, b)
        b .*= -1.0
        α = (1 + γ) * adotAa
        push!(Ellipsoids, EllipsoidBBIS(A, b, α))
    end
    return Ellipsoids
end

##

"""
CRMprod(x₀, Ellipsoids)
Uses cCRM to find a point into intersection of two EllipsoidBBIS 
"""
function CRMprod(x₀::Vector,
    Ellipsoids::Vector{EllipsoidBBIS};
    kwargs...)
    Projections = Function[]
    for ell in Ellipsoids
        push!(Projections, x -> Proj_Ellipsoid(x, ell))
    end
    return CRMprod(x₀, Projections; kwargs...)
end


"""
CARMprod(x₀, Ellipsoids)
Uses CARMprod to find a point into intersection of p EllipsoidBBIS 
"""
function CARMprod(x₀::Vector,
    Ellipsoids::Vector{EllipsoidBBIS};
    kwargs...)
    Projections = Function[]
    for ell in Ellipsoids
        push!(Projections, x -> ApproxProj_Ellipsoid(x, ell))
    end
    return CRMprod(x₀, Projections; kwargs...)
end



"""
MAP(x₀, Ellipsoids)
Uses MAP to find  a point into intersection of p EllipsoidBBIS 
"""
function MAPprod(x₀::Vector,
    Ellipsoids::Vector{EllipsoidBBIS};
    kwargs...)
    Projections = Function[]
    for ell in Ellipsoids
        push!(Projections, x -> Proj_Ellipsoid(x, ell))
    end
    return MAPprod(x₀, Projections; kwargs...)
end



"""
MAAP(x₀, Ellipsoids)
Uses MAAP to find  a point into intersection of p EllipsoidBBIS 
"""
function MAAPprod(x₀::Vector,
    Ellipsoids::Vector{EllipsoidBBIS};
    kwargs...)
    Projections = Function[]
    for ell in Ellipsoids
        push!(Projections, x -> ApproxProj_Ellipsoid(x, ell))
    end
    return MAPprod(x₀, Projections; kwargs...)
end


##



"""
TestEllipsoids()

"""
function TestEllipsoids(; n::Int = 100,
    num_sets::Int = 10,
    samples::Int = 1,
    ε::Real = 1e-6,
    itmax::Int = 1000,
    restarts::Int = 1,
    print_file::Bool = false,
    method::Vector{Symbol} = [:CARMprod, :CRMprod],
    bench_time::Bool = false)
    # X = R^n
    # Defines DataFrame for Results
    dfResults, dfFilenames = createDaframes(method, useapprox)
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
            for mtd in method
                func = eval(mtd)
                filename = savename("Are21", (mtd = mtd, time = timenow), "csv", sort = false)
                print_file ? filedir = datadir("sims", filename) : filedir = ""
                results = func(x₀, Ellipsoids, itmax = itmax, EPSVAL = ε, gap_distance = true, filedir = filedir)
                elapsed_time = 0.0
                if bench_time
                    t = @benchmark $func($x₀, $Ellipsoids, itmax = $itmax, EPSVAL = $ε, gap_distance = true, filedir = $filedir)
                    elapsed_time = (mean(t).time) * 1e-9
                end
                push!(dfrow, results.iter_total)
                push!(dfrow, elapsed_time)
                push!(dfrowFilename, filedir)
                if useapprox
                    Ellipsoids[begin][:useapprox] = true
                    mtd = Symbol("Approx" * String(mtd))
                    filename = savename("AABBIS21", (mtd = mtd, time = timenow), "csv", sort = false)
                    print_file ? filedir = datadir("sims", filename) : filedir = ""
                    results = func(x₀, Ellipsoids, itmax = itmax, EPSVAL = ε, gap_distance = true, filedir = filedir)
                    elapsed_time = 0.0
                    if bench_time
                        t = @benchmark $func($x₀, $Ellipsoids, itmax = $itmax, EPSVAL = $ε, gap_distance = true, filedir = $filedir)
                        elapsed_time = (mean(t).time) * 1e-9
                    end
                    push!(dfrow, results.iter_total)
                    push!(dfrow, elapsed_time)
                    push!(dfrowFilename, filedir)
                    Ellipsoids[begin][:useapprox] = false

                end
            end
            push!(dfResults, dfrow)
            push!(dfFilenames, dfrowFilename)
        end
    end
    return dfResults, dfFilenames
end

##

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







