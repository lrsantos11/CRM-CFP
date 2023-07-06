##################################################################
### Methods
##################################################################


"""
    CRMprod(x₀, Ellipsoids)

Circumcentered-Reflection method for the product space for ellipsoids

"""
function CRMprod(x₀::AbstractVector,
              Ellipsoids::AbstractVector{EllipsoidCRM}; kargs...)
              
    xCRMprod = [x₀ for _ in Ellipsoids]
    ApproxProjections(X) = ApproxProjectEllipsoids_ProdSpace(X, Ellipsoids)
    return CRMprod(xCRMprod, ApproxProjections; kargs...)
end
"""
    DRMprod(x₀, Ellipsoids)

Douglas-Rachford method for the product space for ellipsoids

"""
function DRMprod(x₀::AbstractVector,
    Ellipsoids::AbstractVector{EllipsoidCRM}; kargs...)
    xDRMprod = [x₀ for _ in Ellipsoids]
    xDRMprod[end] = [0.0 for _ in x₀]
    ApproxProjections(X) = ApproxProjectEllipsoids_ProdSpace(X, Ellipsoids)
    return DRMprod(xDRMprod, ApproxProjections; kargs...)
end

ϵ1(k) = 1 / sqrt(k + 1)
ϵ2(k) = 1 / (k + 1)


"""
    PACA(x₀, Ellipsoids)

Perturbed Approximate Circumcenter Algorithm for ellipsoids

"""
function PACA(x₀::AbstractVector,
    Ellipsoids::AbstractVector{EllipsoidCRM},
    ϵ::Function; 
    kargs...)
    FuncEllipsoids = Function[x -> eval_EllipsoidCRM(x, ell) for ell in Ellipsoids]
    ∂Ellipsoids = Function[x -> gradient_EllipsoidCRM(x, ell) for ell in Ellipsoids]
    return PACA(x₀, FuncEllipsoids, ∂Ellipsoids, ϵ; kargs...)
end

PACA1(x₀::AbstractVector,
    Ellipsoids::AbstractVector{EllipsoidCRM}; 
    kargs...) = PACA(x₀, Ellipsoids, ϵ1 ; kargs...)


PACA2(x₀::AbstractVector,
    Ellipsoids::AbstractVector{EllipsoidCRM};
    kargs...) = PACA(x₀, Ellipsoids, ϵ2; kargs...)


"""
    MCSPM(x₀, Ellipsoids)

Modified Cyclic Subgradient Projection method for ellipsoids from [Pierro:1988]

[Pierro:1988] A. R. De Pierro and A. N. Iusem, “A finitely convergent ‘row-action’ method for the convex feasibility problem,” Appl Math Optim, vol. 17, no. 1, Art. no. 1, Jan. 1988, doi: 10.1007/BF01448368.

"""
function MCSPM(x₀::AbstractVector,
    Ellipsoids::AbstractVector{EllipsoidCRM},
    ϵ::Function; 
    kargs...)
    FuncEllipsoids = Function[x -> eval_EllipsoidCRM(x, ell) for ell in Ellipsoids]
    ∂Ellipsoids = Function[x -> gradient_EllipsoidCRM(x, ell) for ell in Ellipsoids]
    return MCSPM(x₀, FuncEllipsoids, ∂Ellipsoids, ϵ; kargs...)
end


MCSPM1(x₀::AbstractVector,
    Ellipsoids::AbstractVector{EllipsoidCRM};
    kargs...) = MCSPM(x₀, Ellipsoids, ϵ1; kargs...)


MCSPM2(x₀::AbstractVector,
    Ellipsoids::AbstractVector{EllipsoidCRM};
    kargs...) = MCSPM(x₀, Ellipsoids, ϵ2; kargs...)



"""
    MSSPM(x₀, Ellipsoids)

Modified simultaneous Subgradient Projection method for ellipsoids from [Iusem:1986]

[Iusem:1986] A. N. Iusem and L. Moledo, “A finitely convergent method of simultaneous subgradient projections for the convex feasibility problem,” Matemática Aplicada e Computacional, vol. 5, no. 2, pp. 169–184, 1986.

"""
function MSSPM(x₀::AbstractVector,
    Ellipsoids::AbstractVector{EllipsoidCRM},
    ϵ::Function;
    kargs...)
    FuncEllipsoids = Function[x -> eval_EllipsoidCRM(x, ell) for ell in Ellipsoids]
    ∂Ellipsoids = Function[x -> gradient_EllipsoidCRM(x, ell) for ell in Ellipsoids]
    return MSSPM(x₀, FuncEllipsoids, ∂Ellipsoids, ϵ; kargs...)
end


MSSPM1(x₀::AbstractVector,
    Ellipsoids::AbstractVector{EllipsoidCRM};
    kargs...) = MSSPM(x₀, Ellipsoids, ϵ1; kargs...)


MSSPM2(x₀::AbstractVector,
    Ellipsoids::AbstractVector{EllipsoidCRM};
    kargs...) = MSSPM(x₀, Ellipsoids, ϵ2; kargs...)




##
# Prototype of the function TestEllipsoidsRn
Random.seed!(1234)
itmax = 2_000
n = 100
m = 50
p = 2 * inv(n)
λ = 1.15
ε = 1e-12
Ellipsoids, _ = SampleEllipsoids(n, m, p, λ = λ)
x₀ = InitialPoint_EllipsoidCRM(Ellipsoids, n)
ϵ(k) =  1/sqrt(k+1)
resultsPACA1 = PACA1(x₀, Ellipsoids, EPSVAL=ε, itmax=itmax)
resultsMCSPM1 = MCSPM1(x₀, Ellipsoids, EPSVAL=ε, itmax=itmax)
resultsMSSPM1 = MSSPM1(x₀, Ellipsoids, EPSVAL=ε, itmax=itmax)

ϵ(k) = 1 / (k + 1)
resultsPACA2 = PACA2(x₀, Ellipsoids, EPSVAL=ε, itmax=itmax)
resultsMCSPM2 = MCSPM2(x₀, Ellipsoids, EPSVAL=ε, itmax=itmax)
resultsMSSPM2 = MSSPM2(x₀, Ellipsoids, EPSVAL=ε, itmax=itmax)
resultsCRMprod = CRMprod(x₀, Ellipsoids, EPSVAL=ε, itmax=itmax)
resultsDRMprod = DRMprod(x₀, Ellipsoids, EPSVAL=ε, itmax=itmax)




@info "PACA ε(k) = 1/√(k+1)" resultsPACA1.iter_total 
@info "PACA ε(k) = 1/(k+1)" resultsPACA2.iter_total
@info "MCSPM ε(k) = 1/√(k+1)" resultsMCSPM1.iter_total
@info "MCSPM ε(k) = 1/√(k+1)" resultsMCSPM2.iter_total
@info "MSSPM ε(k) = 1/√(k+1)" resultsMSSPM1.iter_total
@info "MSSPM ε(k) = 1/√(k+1)" resultsMSSPM2.iter_total
@info "CRMprod" resultsCRMprod.iter_total
@info "DRMprod" resultsDRMprod.iter_total


##
"""
TestEllipsoids()

"""
function TestEllipsoidsRn(
    n::Int, # dimension
    m::Int; # number of ellipsoids
    samples::Int=1,
    λ::Real=1.15,
    ε::Real=1e-6,
    itmax::Int=200_000,
    restarts::Int=1,
    print_file::Bool=false,
    Methods::AbstractVector{Symbol}=[:PACA1, :PACA2, :DRMprod, :CRMprod, :MCSPM1,  :MCSPM2, :MSSPM1, :MSSPM2],
    bench_time::Bool=false,
    gap_distance=false,
    verbose::Bool=false)
    # X = R^n
    # Defines DataFrame for Results
    dfResults, dfFilenames = createDataframes(Methods, projections=true)
    # Fix Random
    Random.seed!(123)
    # Sparsity of first ellipsoid
    p = 2 * inv(n)
    for j = 1:samples
        ℰ, _ = SampleEllipsoids(n, m, p, λ=λ)
        for i = 1:restarts
            x₀ = InitialPoint_EllipsoidCRM(ℰ, n)
            prob_name = savename((Prob=j, Rest=i, dim=n, numsets=m, lambda=λ); equals="", sort=false)
            @info prob_name
            timenow = Dates.now()
            dfrow = []
            dfrowFilename = []
            push!(dfrow, prob_name)
            push!(dfrowFilename, prob_name)
            for mtd in Methods
                func = eval(mtd)
                @info "method: $(string(func))"
                filename = savename("BBILS24", (time=timenow, mtd=mtd, prob=prob_name,), "csv", sort=false)
                print_file ? filedir = datadir("sims", filename) : filedir = ""
                results = func(x₀, ℰ, EPSVAL=ε, filedir=filedir, itmax=itmax, verbose=verbose)
                iszero(results.proj_total) ? results.proj_total = div(results.iter_total,2)*m : nothing
                verbose && @info mtd, results.iter_total, results.proj_total, results.final_tol
                elapsed_time = 0.0
                if bench_time
                    t = @benchmark $func($x₀, $ℰ, EPSVAL=$ε, filedir=$filedir, itmax=$itmax)
                    elapsed_time = (mean(t).time) * 1e-9
                end
                push!(dfrow, results.iter_total)
                push!(dfrow, results.proj_total) 
                push!(dfrow, results.final_tol)
                push!(dfrow, elapsed_time)
                push!(dfrowFilename, filedir)
            end
            push!(dfResults, dfrow)
            push!(dfFilenames, dfrowFilename)
        end
    end
    return dfResults, dfFilenames
end
##
#################################################################
### Run tests of Paper
##################################################################
Methods = Symbol[:PACA1, :PACA2, :CRMprod, :MCSPM1, :MCSPM2, :MSSPM1, :MSSPM2]

function run_Tests(Methods; write_results::Bool=true, samples=20)
    dimensions = [
        (20, 5)
        (20, 10)
        (20, 20)
        (50, 5)
        (50, 10)
        (50, 20)
        (100, 5)
        (100, 10)
        (100, 20)
    ]
    ε = 1e-6
    λ = 1.15
    itmax = 300_000
    bench_time = true
    dfResultsFinal, _ = createDataframes(Methods, projections=true)
    dfResultsFinal[!, :dim] = Int[]
    dfResultsFinal[!, :num_sets] = Int[]
    ##
    for (dim, num_sets) ∈ dimensions
        @info "dim: $dim, num_sets: $num_sets"
        dfResults, _ = TestEllipsoidsRn(dim, num_sets, Methods=Methods, bench_time=bench_time, verbose=false, itmax=itmax, samples=samples, λ=λ)
        dfResults[!, :dim] .= dim
        dfResults[!, :num_sets] .= num_sets
        append!(dfResultsFinal, dfResults)
    end
    ### Write results
    if write_results    
        timenow = Dates.now()
        file_name = savename("BBILS24_EllipsoidsCFP", (@dict timenow), "csv")
        CSV.write(datadir("sims", file_name), dfResultsFinal)
    end
end
run_Tests(Methods)
##

##
#################################################################
### Make Performance Profiles.
##################################################################
# file_name = "BBILS24_EllipsoidsCFP_timenow=2023-06-22T01:00:21.252.csv"
# print_perprof = false
# if print_perprof
#     include(srcdir("Plots_util.jl"))
#     pgfplotsx()
#     dfResultsPerprof = CSV.read(datadir("sims", file_name), DataFrame)

#     # Total  number of projections
#     T_Projs = Matrix{Float64}(dfResultsPerprof[:, [:SucCentCRM_Cyclic_projs, :SucCentCRM_AlmostViolatedFunc_projs, :SePM_projs, :CRMprod_projs]])
#     T_Projs[findall(row -> row >= itmax, T_Projs)] .= Inf
#     T_Projs = T_Projs[61:end, :] #remove first 60 points

#     perprof1 = performance_profile(PlotsBackend(),
#         T_Projs,
#         # logscale = false,
#         ["Alg1", "Alg3", "SePM", "CRMprod"],
#         legend=:bottomright, framestyle=:box,
#         linestyles=[:dash, :solid, :dot, :dashdot],
#     )
#     ylabel!("Fraction of problems solved")
#     ticks = [0, 2, 4, 6, 8, 10, 12]
#     xticks!(perprof1, ticks, [L"2^{%$(i)}" for i in ticks])
#     title!("Performance Profile -- Total projections comparison")
#     perprof1_file_name = "BBILS23_Ellipsoids_Perprof_TotalProjections.pdf"
#     savefig(perprof1, plotsdir(perprof1_file_name))
#     perprof1

#     ##
#     # Total CPU time

#     T_CPU = Matrix{Float64}(dfResultsPerprof[:, [:SucCentCRM_Cyclic_elapsed, :SucCentCRM_AlmostViolatedFunc_elapsed, :SePM_elapsed, :CRMprod_elapsed]])
#     T_CPU[findall(row -> row >= itmax, T_CPU)] .= Inf
#     T_CPU = T_CPU[61:end, :] #remove first 60 points

#     perprof2 = performance_profile(PlotsBackend(),
#         T_CPU,
#         # logscale = false,
#         ["Alg1", "Alg3", "SePM", "CRMprod"],
#         legend=:bottomright, framestyle=:box,
#         linestyles=[:dash, :solid, :dot, :dashdot],
#     )
#     ylabel!("Fraction of problems solved")
#     ticks = [0, 2, 4, 6, 8, 10, 12]
#     xticks!(perprof2, ticks, [L"2^{%$(i)}" for i in ticks])
#     title!("Performance Profile -- CPU time comparison")
#     perprof2_file_name = "BBILS23_Ellipsoids_Perprof_CPUTime.pdf"
#     savefig(perprof2, plotsdir(perprof2_file_name))
#     perprof2
#     for file in [perprof1_file_name, perprof2_file_name]
#         cp(plotsdir(file),  "../../../Draft/New/Successive-cCRM/"*file, force=true)
#     end
    

# end


##
#################################################################
### Make Tables
##################################################################

# print_tables = true

# if print_tables
#     Methods = [:PACA1, :PACA2, :CRMprod,  :MSSPM1, :MSSPM2]
#     dfResultsTables = CSV.read(datadir("sims", file_name), DataFrame)
#     for type  in ["_projs","_elapsed"]
#         @info("Results for $(type[2:end])")
#         @info describe(dfResultsTables[:, string.(Methods).*type], :mean, :std, :median, :min, :max)
#         df_tabela = DataFrame(n=Int[], m=Int[])
#         for method in Methods
#             df_tabela[!, string(method)*"_mean"] = Float64[]
#             df_tabela[!, string(method)*"_std"] = Float64[]
#         end

#         gdf = groupby(dfResultsTables, [:dim, :num_sets])
#         for k in keys(gdf)
#             df = describe(gdf[k][:,string.(Methods).*type], :mean, :std)    
#             df_row = [] 
#             push!(df_row, k[1], k[2])
#             for row in eachrow(df)
#                 push!(df_row, row[2], row[3])
#             end 
#             push!(df_tabela,df_row)
#         end
#         @info("Total $(type[2:end])")
#         @info df_tabela
#     end
# end

##

# using BenchmarkTools
# vₖ = [copy(x₀) for i in 1:m]
# xPACA = copy(x₀)
# Functions =     Function[x -> eval_EllipsoidCRM(x, ell) for ell in Ellipsoids]
# Subgrads = Function[x -> gradient_EllipsoidCRM(x, ell) for ell in Ellipsoids]

# ϵₖ = 1e-4
 
# @btime for i in 1:m
#     copyto!($vₖ[i], computevₖⁱ($xPACA, $Functions[i], $Subgrads[i], ϵ=ϵₖ))
# end

# @btime copyto!(vₖ, [computevₖⁱ(xPACA, Functions[i], Subgrads[i], ϵ=ϵₖ) for i in 1:m])

# @btime Threads.@threads for i in 1:m
#     copyto!($vₖ[i], computevₖⁱ($xPACA, $Functions[i], $Subgrads[i], ϵ=ϵₖ))
# end

# @btime computevₖ!($vₖ, $xPACA, $Functions, $Subgrads, m,ϵ=ϵₖ);


