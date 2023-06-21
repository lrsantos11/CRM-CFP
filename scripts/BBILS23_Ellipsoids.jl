##################################################################
### Methods
##################################################################


"""
    SucCentCRM_Cyclic(x₀, Ellipsoids)

Successive Centralized Cirumcentered-Reflection method using Cyclic control sequence
"""
function SucCentCRM_Cyclic(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM}; kargs...)
    Projections = Function[x -> Proj_Ellipsoid(x, ell) for ell in Ellipsoids]
    return SucCentCRM_Cyclic(x₀, Projections; kargs...)
end

"""
    SucCentCRM_AlmostViolatedDist(x₀, Ellipsoids)
"""
function SucCentCRM_AlmostViolatedDist(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM}; kargs...)
    Projections = Function[x -> Proj_Ellipsoid(x, ell) for ell in Ellipsoids]
    return SucCentCRM_AlmostViolatedDist(x₀, Projections; kargs...)
end

"""
    SucCentCRM_AlmostViolatedFunc(x₀, Ellipsoids)
"""
function SucCentCRM_AlmostViolatedFunc(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM}; kargs...)
    Projections = Function[x -> Proj_Ellipsoid(x, ell) for ell in Ellipsoids]
    FunctionEval = Function[x -> eval_EllipsoidCRM(x, ell) for ell in Ellipsoids]
    return SucCentCRM_AlmostViolatedFunc(x₀, Projections, FunctionEval; kargs...)
end


"""
    SePM(x₀, Ellipsoids)
"""
function SePM(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM}; kargs...)
    Projections = Function[x -> Proj_Ellipsoid(x, ell) for ell in Ellipsoids]
    return SePM(x₀, Projections; kargs...)
end

"""
    CRMProd(x₀, Ellipsoids)
"""
function CRMprod(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM}; kargs...)
    m = length(Ellipsoids)
    Projections = Function[x -> Proj_Ellipsoid(x, ell) for ell in Ellipsoids]
    Results = CRMprod(x₀, Projections; kargs...)
    Results.iter_total = div(Results.iter_total, 2)
    Results.proj_total = m * Results.iter_total
    return Results
end



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
    Methods::Vector{Symbol}=[:SucCentCRM_Cyclic, :SucCentCRM_AlmostViolatedDist, :SucCentCRM_AlmostViolatedFunc, :SePM, :CRMprod],
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
                filename = savename("BBILS23", (time=timenow, mtd=mtd, prob=prob_name,), "csv", sort=false)
                print_file ? filedir = datadir("sims", filename) : filedir = ""
                results = func(x₀, ℰ, EPSVAL=ε, filedir=filedir, itmax=itmax, verbose=verbose)
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
samples = 20
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
Methods = [:SucCentCRM_Cyclic, :SucCentCRM_AlmostViolatedDist, :SucCentCRM_AlmostViolatedFunc, :SePM, :CRMprod]
dfResultsFinal, _ = createDataframes(Methods, projections=true)
dfResultsFinal[!, :dim] = Int[]
dfResultsFinal[!, :num_sets] = Int[]
##
for (dim, num_sets) ∈ dimensions
    dfResults, _ = TestEllipsoidsRn(dim, num_sets, Methods=Methods, bench_time=bench_time, verbose=false, itmax=itmax, samples=samples, λ=λ)
    dfResults[!, :dim] .= dim
    dfResults[!, :num_sets] .= num_sets
    append!(dfResultsFinal, dfResults)
end

##
#################################################################
### Write results
##################################################################
timenow = Dates.now()
file_name = savename("BBILS23_EllipsoidsCFP", (@dict timenow), "csv")
CSV.write(datadir("sims", file_name), dfResultsFinal)


##
#################################################################
### Make Performance Profiles.
##################################################################
print_perprof = false
if print_perprof
    include(srcdir("Plots_util.jl"))
    pgfplotsx()
    file_name = "BBILS23_EllipsoidsCFP_timenow=2023-05-10T07:03:09.048.csv"
    dfResultsPerprof = CSV.read(datadir("sims", file_name), DataFrame)

    # Total  number of projections
    T_Projs = Matrix{Float64}(dfResultsPerprof[:, [:SucCentCRM_Cyclic_projs, :SucCentCRM_AlmostViolatedFunc_projs, :SePM_projs, :CRMprod_projs]])
    T_Projs[findall(row -> row >= itmax, T_Projs)] .= Inf
    T_Projs = T_Projs[61:end, :] #remove first 60 points

    perprof1 = performance_profile(PlotsBackend(),
        T_Projs,
        # logscale = false,
        ["Alg1", "Alg3", "SePM", "CRMprod"],
        legend=:bottomright, framestyle=:box,
        linestyles=[:dash, :solid, :dot, :dashdot],
    )
    ylabel!("Fraction of problems solved")
    ticks = [0, 2, 4, 6, 8, 10, 12]
    xticks!(perprof1, ticks, [L"2^{%$(i)}" for i in ticks])
    title!("Performance Profile -- Total projections comparison")
    perprof1_file_name = "BBILS23_Ellipsoids_Perprof_TotalProjections.pdf"
    savefig(perprof1, plotsdir(perprof1_file_name))
    perprof1

    ##
    # Total CPU time

    T_CPU = Matrix{Float64}(dfResultsPerprof[:, [:SucCentCRM_Cyclic_elapsed, :SucCentCRM_AlmostViolatedFunc_elapsed, :SePM_elapsed, :CRMprod_elapsed]])
    T_CPU[findall(row -> row >= itmax, T_CPU)] .= Inf
    T_CPU = T_CPU[61:end, :] #remove first 60 points

    perprof2 = performance_profile(PlotsBackend(),
        T_CPU,
        # logscale = false,
        ["Alg1", "Alg3", "SePM", "CRMprod"],
        legend=:bottomright, framestyle=:box,
        linestyles=[:dash, :solid, :dot, :dashdot],
    )
    ylabel!("Fraction of problems solved")
    ticks = [0, 2, 4, 6, 8, 10, 12]
    xticks!(perprof2, ticks, [L"2^{%$(i)}" for i in ticks])
    title!("Performance Profile -- CPU time comparison")
    perprof2_file_name = "BBILS23_Ellipsoids_Perprof_CPUTime.pdf"
    savefig(perprof2, plotsdir(perprof2_file_name))
    perprof2
    for file in [perprof1_file_name, perprof2_file_name]
        cp(plotsdir(file),  "../../../Draft/New/Successive-cCRM/"*file, force=true)
    end
    

end


##
#################################################################
### Make Tables
##################################################################

print_tables = false

if print_tables

    file_name = "BBILS23_EllipsoidsCFP_timenow=2023-05-10T07:03:09.048.csv"
    dfResultsTables = CSV.read(datadir("sims", file_name), DataFrame)
    for tipo  in ["_projs","_elapsed"]
        @show describe(dfResultsTables[:,string.(Methods).*tipo], :mean, :std, :median, :min, :max)
        df_tabela = DataFrame(n=Int[],m=Int[], Alg1_mean=Float64[], Alg1_std=Float64[], Alg_2_mean = Float64[],  Alg_2_std = Float64[], Alg3_mean = Float64[], Alg3_std = Float64[], SePM_mean = Float64[],  SePM_std = Float64[], CRMprod_mean = Float64[],  CRMprod_std = Float64[])

        gdf = groupby(dfResultsTables, [:dim, :num_sets])
        for k in keys(gdf)
            df = describe(gdf[k][:,string.(Methods).*tipo], :mean, :std)    
            df_row = [] 
            push!(df_row, k[1], k[2])
            for row in eachrow(df)
                push!(df_row, row[2], row[3])
            end 
            push!(df_tabela,df_row)
        end
        println("Total number of $tipo")
        @show df_tabela
    end
end