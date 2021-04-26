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
function TestPolyhedral(;ninit::Int64 = 200,samples::Int64 = 1,
        ε::Float64 = 1e-5,itmax::Int64 = 1000, restarts = 1,print_file::Bool=false)
    # X = R^n
    # Defines DataFrame for Results
    dfResults= DataFrame(Problem=String[],CRMit=Int[],DRMit=Int[],MAPit=Int[])
    dfFilenames= DataFrame(Problem=String[],CRMfilename=String[],DRMfilename=String[],MAPfilename=String[])

    # Fix Random
    Random.seed!(1)
    for j in 1:samples
        # n = rand(ninit:2*ninit)
        n = ninit
        # Generates Matrix m × n  with m < n
        m = rand(1:n-1)
        A  = randn(m,n)
        m, n = size(A)
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
            x₀ = StartingPoint(n)
            prob_name  = String("Problem$j"*"Restart$i")
            timenow= Dates.now()
            dfrow = []
            dfrowFilename = []
            push!(dfrow,prob_name)
            push!(dfrowFilename,prob_name)
            methods = [:CRMprod, :DRMprod, :MAPprod]
            for mtd in methods
                func = eval(mtd) 
                filename = savename("BBS20Sec42",(mtd=mtd,time=timenow),"csv",sort=false)
                print_file ? filedir = datadir("sims",filename) : filedir = ""
                results  = func(x₀,halfspacesInd,itmax=itmax,EPSVAL=ε,gap_distance=true,filedir=filedir)
                push!(dfrow,results.iter_total)
                push!(dfrowFilename,filedir)
            end    
            push!(dfResults,dfrow)
            push!(dfFilenames,dfrowFilename)
        end
    end
    return dfResults,dfFilenames
end

# Constructs BBS20Fig4_TestPolyhedral
ninit = 200
samples = 1
restarts = 1
ε = 1e-6
itmax = 20000
print_file = true

dfResultsPoly,dfrowFilename = TestPolyhedral(ninit = ninit, samples = samples,itmax=itmax, ε = ε, restarts = restarts,print_file=print_file)

xCRM = readdlm(dfrowFilename.CRMfilename[1])
xDRM = readdlm(dfrowFilename.DRMfilename[1])
xMAP = readdlm(dfrowFilename.MAPfilename[1])
plt_poly = plot(xCRM[1:end,1].+1,xCRM[1:end,2],scale=:log10, label="CRM-prod",
            title="Comparison using Product Space reformulation",
            framestyle = :box,
           xlabel = "Number of iterations (log scale)",
           ylabel = "Gap error (log scale)",
           minorticks=true)
plot!(xDRM[1:end,1].+1,xDRM[1:end,2] .+ 1e-6,scale=:log10, label="DRM-prod",linestyle=:dash,minorticks=false)
plot!(xMAP[1:end,1].+1,xMAP[1:end,2],scale=:log10, label="MAP-prod",linestyle=:dot)
##
savefig(plt_poly,plotsdir("BBS20Fig4_TestPolyhedral.pdf"))
##
# Constructs BBS20Sec4 Table
ninit = 200
samples = 100
restarts = 10
ε = 1e-6
itmax = 20000
print_file = false

dfResultsPoly,dfrowFilename = TestPolyhedral(ninit = ninit, samples = samples,itmax=itmax, ε = ε, restarts = restarts,print_file=print_file)
describe(dfResultsPoly)




