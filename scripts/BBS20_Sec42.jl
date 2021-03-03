"""
This script builds the results and plots presented in Section 4.2 of [^Behling2020]

[^Behling2020] R. Behling, J.-Y. Bello-Cruz, e L.-R. Santos, “On the Circumcentered-Reflection Method 
for the Convex Feasibility Problem”, Numer. Algorithms, jul. 2020, doi: [10.1007/s11075-020-00941-6](https://doi.org/10.1007/s11075-020-00941-6). 
"""

include(scriptsdir("plots_util.jl"))

"""
TestPolyhedral(;n::Int64 = 100,samples::Int64 = 2,
        EPS_VAL::Float64 = 1e-5,printfile::Bool=false,itmax::Int64 = 200)

"""
function TestPolyhedral(;n::Int64 = 200,samples::Int64 = 1,
        EPSVAL::Float64 = 1e-5,itmax::Int64 = 1000, restarts = 1,print_file::Bool=false)
    # X = R^n
    # Fix Random
    Random.seed!(1)
    # Defines DataFrame for Results
    dfResults= DataFrame(Problem=String[],CRMit=Int[],DRMit=Int[],MAPit=Int[])
    for j in 1:samples
        n = rand(n:2*n)
        # Generates Matrix m × n  with m < n
        m = rand(1:n-1)
        rand()
        A  = randn(m,n)
        @show m, n = size(A)
        xbar = StartingPoint(n)
        bbar = A*xbar
        # Number of halfspaces
        non_slater = rand(1:m)
        indices = sort(rand(1:m,non_slater))
        r = rand(non_slater)
        bbar[indices] +=  norm(bbar[indices])*r
        halfspacesInd = ProximableFunction[]
        for index  in eachindex(bbar)
            push!(halfspacesInd,IndHalfspace(A[index,:],bbar[index]))
        end
        for i = 1:restarts
            xzero = StartingPoint(n)
            prob_name  = String("Problem$j"*"Restart$i")
            timenow= Dates.now()
            dfrow = []
            push!(dfrow,prob_name)
            methods = [:CRMprod]#, :DRMprod, :MAPprod]
            for mtd in methods
                filename = savename("BBS20Sec42",(mtd=:CRMprod,time=timenow),"csv",sort=false)
                print_file ? filedir = datadir("sims",filename) : filedir = ""
                resultCRM  = CRMprod(xzero,halfspacesInd,itmax=itmax,EPSVAL=EPSVAL,gap_distance=true,filedir=filedir)
                # @show resultCRM
                push!(dfResults,(prob_name,resultCRM.iter_total,resultCRM.iter_total,resultCRM.iter_total))
            end    
        end
    end
    return dfResults
end
n = 200
samples = 1
restarts = 1
ε = 1e-5
itmax = 20000
mtd=:CRMProd
print_file = true

dfResultsPoly = TestPolyhedral(n = n, samples = samples,itmax=itmax, EPSVAL = ε, restarts = restarts,print_file=print_file)
describe(dfResultsPoly)
##
plt_poly = plot(xCRM[:,1],xCRM[:,2],scale=:log10, label="CRM-prod",
           title="Comparison using Product Space reformulation",
           framestyle = :box,
           xlabel = "Number of iterations (log scale)",
           ylabel = "Gap error (log scale)",
           minorticks=true);
plot!(xDRM[:,1],xDRM[:,2],scale=:log10, label="DRM-prod",linestyle=:dash);
plot!(xMAP[:,1],xMAP[:,2],scale=:log10, label="MAP-prod",linestyle=:dot); # linestyles=[:solid, :dash, :dot])

savefig(plt_poly,plotsdir("BBS20Fig4_TestPolyhedral.pdf"))
@show describe(dfResultsPoly)
plt_poly

# p1 = plot(xCRM[:,1],xCRM[:,2],scale=:log10, label="CRM-prod",
#            title="Comparison using Product Space reformulation",
#            framestyle = :box,
#            xlabel = "Number of iterations (log scale)",
#            ylabel = "Gap error (log scale)",
#            minorticks=true);
# p1 = plot!(xDRM[:,1],xDRM[:,2],scale=:log10, label="DRM-prod",linestyle=:dash);
# p1 = plot!(xMAP[:,1],xMAP[:,2],scale=:log10, label="MAP-prod",linestyle=:dot); # linestyles=[:solid, :dash, :dot])
# # plot(p1)
# # performance_profile(hcat(resultsCRM[:,2],resultsDRM[:,2],resultsMAP[:,2]), ["CRM", "DRM", "MAP"], title="Polyhedral Intersection")
# # savefig(p1, "tables/Polyhedral"*string(now())*".pdf")
# savefig(p1, "../../Draft/New/Two-Any-Convex-Sets/figures/ComparisonCRMProd_DRM-prod_MAP-prod-bw.pdf")
# df = DataFrame()
# df.CRM = resultsCRM[:,2]
# df.DRM = resultsDRM[:,2]
# df.MAP = resultsMAP[:,2]
# describe(df)




