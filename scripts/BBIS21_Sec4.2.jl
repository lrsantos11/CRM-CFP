
"""
This script builds the results and plots presented in Section 4.1 of [^Behling2022]

[^Behling2022] R. Behling, J.-Y. Bello-Cruz, A. Iusem and L.-R. Santos, 
"""
##################################################################
## Basic Functions for and Polyhedral and SOC tests
##################################################################
using BenchmarkTools


function TestPolyhedralSOC(;
    n::Int64=200,
    samples::Int64=10,
    restarts::Int64=1,
    EPSVAL::Float64=1e-6,
    itmax::Int64=2_000,
    t::Float64=0.0)
    # Fix Random
    Random.seed!(10)
    # Defines DataFrame for Results
    dfResults = DataFrame(Problem=String[], CRMprodproj=Int[], cCRMproj=Int[], DRMproj=Int[], MAPproj=Int[])
    #md Defines Indicator Function of Second Order Cone from `ProximalOperators.jl`
    for j in 1:samples
        #########################
        w = StartingPoint(n - 1)
        xSol = [norm(w); w]
        vperp = [-norm(w); w]
        num_rows = rand(div(n, 3):n)
        Apolar, bpolar = generate_polar(vperp, num_rows, xSol - t * vperp)
        Poly = IndPolyhedral(Apolar, bpolar)
        SOC = IndSOC()
        ProjectA(x) = ProjectIndicator(SOC, x)
        ProjectB(x) = ProjectIndicator(Poly, x)
        # Restarts
        for i = 1:restarts
            xzero = StartingPoint(n)
            xzero = ProjectB(xzero)
            # @show SOC(xzero)
            while SOC(xzero) != Inf
                xzero = StartingPoint(n)
                xzero = ProjectB(xzero)
            end
            prob_name = String("Problem$j" * "Restart$i")
            # println(prob_name)
            resultCRM = CRMprod(xzero, ProjectA, ProjectB, itmax=itmax, EPSVAL=EPSVAL, gap_distance=true)
            resultcCRM = centralizedCRM(xzero, ProjectA, ProjectB, itmax=itmax, EPSVAL=EPSVAL, gap_distance=true)
            resultDRM = DRM(xzero, ProjectA, ProjectB, itmax=itmax, EPSVAL=EPSVAL, gap_distance=true)
            resultMAP = MAP(xzero, ProjectA, ProjectB, itmax=itmax, EPSVAL=EPSVAL, gap_distance=true)
            push!(dfResults, (prob_name, resultCRM.iter_total, resultcCRM.iter_total, resultDRM.iter_total, resultMAP.iter_total))
        end
    end
    return dfResults
end

##

function normalized_dist_proximity(xk, ProjectA, ProjectB, xzero)
    norm2xkPA = norm(ProjectA(xk) - xk)^2
    norm2xkPB = norm(ProjectB(xk) - xk)^2
    norm2xzeroPA = norm(ProjectA(xzero) - xzero)^2
    norm2xzeroPB = norm(ProjectB(xzero) - xzero)^2
    return 10 * log10((norm2xkPA + norm2xkPB) / (norm2xzeroPA + norm2xzeroPB))
end




function generate_polar(vec::AbstractArray, num_edges::Int, bound_point::AbstractArray; lbound::Float64=0.1)
    size_vec = size(vec)
    vec_unit = vec ./ norm(vec)
    polar = []
    k = 0
    while k < num_edges
        u_polar = randn(size_vec)
        u_polar /= norm(u_polar)
        if lbound <= dot(u_polar, -vec_unit)
            push!(polar, u_polar)
            k += 1
        end
    end
    Apolar = permutedims(hcat(polar...))
    bpolar = Apolar * bound_point
    return Apolar, bpolar
end

# ##
n = 200
samples = 20
restarts = 5
EPSVAL = 1e-6
itmax = 300
dfResults = TestPolyhedralSOC(n=n, samples=samples, itmax=itmax, EPSVAL=EPSVAL, restarts=restarts)
# perprof = performance_profile(hcat(dfResults.CRMit, dfResults.DRMit, dfResults.MAPit), ["CRM", "DRM", "MAP"],
#     title=L"Performance Profile -- Gap error -- $\varepsilon = 10^{-6}$",
#     legend=:bottomright, framestyle=:box, linestyles=[:solid, :dash, :dot])
# ylabel!("Percentage of problems solved")
# savefig(perprof, plotsdir("BBS20Fig3_PolyhedralSOC.pdf"))
# @show describe(dfResults)
# perprof

##


# n = 200
# w = randn(n - 1)
# xSol = [norm(w); w]
# vperp = [-norm(w); w]
# t = 1
# Apolar, bpolar = generate_polar(vperp, 100, xSol - t*vperp,lbound = 0.1)
# Poly = IndPolyhedral(Apolar, bpolar)
# SOC = IndSOC()
# ProjectA(x) = ProjectIndicator(SOC, x)
# ProjectB(x) = ProjectIndicator(Poly, x)
# xzero = StartingPoint(n)
# ##
# resultcCRM = centralizedCRM(xzero, ProjectA, ProjectB, itmax=itmax, EPSVAL=EPSVAL, gap_distance = true)
# ##
# resultMAP = MAP(xzero, ProjectA, ProjectB, itmax=itmax, EPSVAL=1e-14, gap_distance = true)
