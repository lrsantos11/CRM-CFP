
"""
This script builds the results and plots presented in Section 4.1 of [^Behling2022]

[^Behling2022] R. Behling, J.-Y. Bello-Cruz, A. Iusem and L.-R. Santos, 
"""
##################################################################
## Basic Functions for and Polyhedral and SOC tests
##################################################################


function TestPolyhedralSOC(;
    n::Int64=200,
    samples::Int64=10,
    restarts::Int64=1,
    EPSVAL::Float64=1e-6,
    itmax::Int64=2_000,
    τ::Float64=0.0)
    # Fix Random
    Random.seed!(10)
    # Defines DataFrame for Results
    dfResults = DataFrame(Problem=String[], CRMprod=Int[], cCRM=Int[], DRM=Int[], MAP=Int[])
    #md Defines Indicator Function of Second Order Cone from `ProximalOperators.jl`
    for j in 1:samples
        #########################
        w = StartingPoint(n - 1)
        xSol = [norm(w); w]
        vperp = [-norm(w); w]
        @show num_rows = rand(div(n, 3):n)
        Apolar, bpolar = generate_polar(vperp, num_rows, xSol - τ * vperp)
        Poly = IndPolyhedral(Apolar, bpolar)
        SOC = IndSOC()
        ProjectA(x) = ProjectIndicator(SOC, x)
        ProjectB(x) = ProjectIndicator(Poly, x)
        # Restarts
        for i = 1:restarts
            xzero = StartingPoint(n)
            # @show SOC(xzero)
            while SOC(xzero) != Inf
                xzero = StartingPoint(n)
            end
            prob_name = String("Problem$j" * "Restart$i")
            println(prob_name)
            println("CRMProd")
            resultCRMprod = CRMprod(xzero, ProjectA, ProjectB, itmax=itmax, EPSVAL=EPSVAL, gap_distance=true)
            println("cCRM")
            resultcCRM = centralizedCRM(xzero, ProjectA, ProjectB, itmax=itmax, EPSVAL=EPSVAL, gap_distance=true)
            println("DRM")
            resultDRM = DRM(xzero, ProjectA, ProjectB, itmax=itmax, EPSVAL=EPSVAL, gap_distance=true)
            println("MAP")
            resultMAP = MAP(xzero, ProjectA, ProjectB, itmax=itmax, EPSVAL=EPSVAL, gap_distance=true)
            push!(dfResults, (prob_name, resultCRMprod.iter_total, resultcCRM.iter_total, resultDRM.iter_total, resultMAP.iter_total))
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

##
n = 200
samples = 50
restarts = 4
EPSVAL = 1e-6
itmax = 10_000
τ = 0
# dfResults = TestPolyhedralSOC(n=n, samples=samples, itmax=itmax, EPSVAL=EPSVAL, restarts=restarts)
# timenow = Dates.now()
# CSV.write(datadir("sims", savename("BBIS21_PolyhedralSOC", (@dict timenow EPSVAL τ), "csv")), dfResults)
##

dfResults = CSV.read(datadir("sims", "BBIS21_PolyhedralSOC_EPSVAL=1e-6_timenow=2022-09-12T17:19:12.824.csv"), DataFrame)

perprof = performance_profile(PlotsBackend(), float.(hcat(dfResults.cCRM, dfResults.MAP, dfResults.CRMprod)), ["cCRM", "MAP", "CRMprod",],
    # title=L"Performance Profile -- Gap error -- $\varepsilon = 10^{-6}$",
    legend=:bottomright, framestyle=:box, linestyles=[:solid, :dash, :dot])
ylabel!("Percentage of problems solved")
xticks!(perprof, 0:0.5:2.5, [L"2^{0}", L"2^{0.5}", L"2^{1}", L"2^{1.5}", L"2^2", L"2^{2.5}"])
savefig(perprof, plotsdir("BBIS21_PolyhedralSOC_Perprof_tau=$(τ).pdf.pdf"))
@show describe(dfResults, :mean, :std, :median, :min, :max)[[3, 5, 2], :]
perprof

##
n = 200
samples = 50
restarts = 4
EPSVAL = 1e-6
itmax = 10_000
τ = .25
# dfResults = TestPolyhedralSOC(n=n, samples=samples, itmax=itmax, EPSVAL=EPSVAL, restarts=restarts, τ = τ)
# timenow = Dates.now()
# CSV.write(datadir("sims", savename("BBIS21_PolyhedralSOC", (@dict timenow EPSVAL τ), "csv")), dfResults)
##

dfResults2 = CSV.read(datadir("sims", "BBIS21_PolyhedralSOC_EPSVAL=1e-6_timenow=2022-09-15T23:00:13.669_τ=0.25.csv"), DataFrame)

perprof2 = performance_profile(PlotsBackend(), float.(hcat(dfResults2.cCRM, dfResults2.MAP, dfResults2.CRMprod)), ["cCRM", "MAP", "CRMprod",],
    # title=L"Performance Profile -- Gap error -- $\varepsilon = 10^{-6}$",
    legend=:bottomright, framestyle=:box, linestyles=[:solid, :dash, :dot])
ylabel!("Percentage of problems solved")
xticks!(perprof2, 0:0.5:2.5, [L"2^{0}", L"2^{0.5}", L"2^{1}", L"2^{1.5}", L"2^2", L"2^{2.5}"])
savefig(perprof2, plotsdir("BBIS21_PolyhedralSOC_Perprof_tau=$(τ).pdf"))
@show describe(dfResults2, :mean, :std, :median, :min, :max)[[3, 5, 2], :]
perprof2




##

# n = 200
# Random.seed!(10)
# w = StartingPoint(n - 1)
# xSol = [norm(w); w]
# vperp = [-norm(w); w]
# num_rows = rand(div(n, 3):n)
# Apolar, bpolar = generate_polar(vperp, num_rows, xSol - t * vperp)
# Poly = IndPolyhedral(Apolar, bpolar)
# SOC = IndSOC()
# ProjectA(x) = ProjectIndicator(SOC, x)
# ProjectB(x) = ProjectIndicator(Poly, x)
# xzero = StartingPoint(n)
# while SOC(xzero) != Inf
#     xzero = StartingPoint(n)
# end
# ##
# resultcCRM = centralizedCRM(xzero, ProjectA, ProjectB, itmax=itmax, EPSVAL=EPSVAL, gap_distance=true)
# # ##
# resultMAP = MAP(xzero, ProjectA, ProjectB, itmax=itmax, EPSVAL=1e-14, gap_distance = true)
