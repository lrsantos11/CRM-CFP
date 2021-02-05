"""
This script builds the results and plots presented in Section 4.1 of [^Behling2020]

[^Behling2020] R. Behling, J.-Y. Bello-Cruz, e L.-R. Santos, “On the Circumcentered-Reflection Method 
for the Convex Feasibility Problem”, Numer. Algorithms, jul. 2020, doi: [10.1007/s11075-020-00941-6](https://doi.org/10.1007/s11075-020-00941-6). 
"""

include(scriptsdir("plots_util.jl"))

function TestAffineSOC(;n::Int64 = 200,samples::Int64 = 10, restarts::Int64 = 1,
        EPSVAL::Float64 = 1e-6,itmax::Int64 = 2_000)
    # Fix Random
    Random.seed!(10)
    # Defines DataFrame for Results
    dfResults= DataFrame(Problem=String[],CRMit=Int[],DRMit=Int[],MAPit=Int[])
    #md Defines Indicator Function of Second Order Cone from `ProximalOperators.jl`
    SOC = IndSOC()
    for j in 1:samples
        #########################
        # Generate Subspace
        # affine = true
        # cmax = rand(1:ceil(Integer,n/10,))
        # CC  = GenerateSamples(n,affine)
        # #  Read Files
        m = rand(1:n-1)
        A =  randn(m,n)
        b = randn(m)
        w  = (A\b)[2:end]
        b = A*[norm(w); w]
        Affine = IndAffine(A,b)
        ProjectA(x) =  ProjectIndicator(SOC,x)
        ProjectB(x) =  ProjectIndicator(Affine,x)
        # Restarts
        for i = 1:restarts
            xzero = StartingPoint(n)
            xzero = ProjectB(xzero)
            # @show SOC(xzero)
            while SOC(xzero) != Inf
                xzero = StartingPoint(n)
                xzero = ProjectB(xzero)
            end
            prob_name  = String("Problem$j"*"Restart$i")
            # println(prob_name)
            resultCRM  = CRM(xzero,ProjectA,ProjectB,itmax=itmax,EPSVAL=EPSVAL,gap_distance=true)
            resultDRM  = DRM(xzero,ProjectA,ProjectB,itmax=itmax,EPSVAL=EPSVAL,gap_distance=true)
            resultMAP  = MAP(xzero,ProjectA,ProjectB,itmax=itmax,EPSVAL=EPSVAL,gap_distance=true)
            push!(dfResults,(prob_name,resultCRM.iter_total,resultDRM.iter_total,resultMAP.iter_total))
        end
    end
    return dfResults
end


n = 200
samples = 100
restarts = 10
ε = 1e-6
itmax = 2000
dfResults = TestAffineSOC(n = n, samples = samples,itmax=itmax, EPSVAL = ε, restarts = restarts)
perprof = performance_profile(hcat(dfResults.CRMit,dfResults.DRMit,dfResults.MAPit), ["CRM", "DRM", "MAP"],
    title=L"Performance Profile -- Gap error -- $\varepsilon = 10^{-6}$",
    legend = :bottomright, framestyle = :box, linestyles=[:solid, :dash, :dot])
ylabel!("Percentage of problems solved")
savefig(perprof,plotsdir("BBS20Fig3_AffineSOC.pdf"))
@show describe(dfResults)
perprof

