using DrWatson
@quickactivate "CRM-CFP"

using Test
# Include the CRM-CFP.jl file
include(srcdir("CRM-CFP.jl"))

function print_results(method, num_rows, num_cols, xSol, num_blocks, results, elap_time, c, df_results)
    dist_sol = norm(xSol - results.xApprox) / norm(xSol)
    @info "Solving system with size ($num_rows,$num_cols) and with coherence $(c)"
    @info "Method: $(method) with $num_blocks blocks"
    @info "Distance to solution:  $(dist_sol)"
    @info "Total Projections:  $(results.iter_total)"
    @info "Wall-clock time for $(method) with $num_blocks blocks: $(elap_time)"
    println("="^80)
    push!(df_results, (method=method, num_blocks=num_blocks, num_rows=num_rows, num_cols=num_cols, coherence=c, distance=dist_sol, total_projections=results.iter_total, wall_clock_time=elap_time))
end


##
@testset "First Example: Random Matrix" begin
    df_results = DataFrame()
    Random.seed!(42)
    T = Float64
    df_results = DataFrame(:method => Symbol[], :num_blocks => Int[], :num_rows => Int[], :num_cols => Int[], :coherence => Float64[], :distance => Float64[], :total_projections => Int[], :wall_clock_time => Float64[])
   
    for num_rows ∈ [5_000, 7_500, 10_000, 12_500],  num_cols ∈ [100, 250, 500], c ∈ T[0.0, 0.1, 0.2]

        # num_rows, num_cols = 20_000, 1_000
        # A = randn(T, num_rows, num_cols)
        A = (1 - c) * randn(T, num_rows, num_cols) + c * ones(T, num_rows, num_cols)
        w = randn(T, num_rows)
        xSol = A' * w
        normalize!(xSol)
        b = A * xSol
        x₀ = zeros(T, num_cols)
        num_blocks = fld(num_rows,  num_cols) + 1
        @info "Solving system with size $(size(A)) and using Julia's QR"
        println("="^80)
        println("="^80)
        for  method in [:pCRM, :CRM] 
            
            Projections = projection_block_QR(A, b, num_blocks)
            func = eval(method)
            results = func(x₀, Projections, num_blocks, itmax=100_000, EPSVAL=1e-5)
            elap_time = @belapsed $(func)($(x₀), $(Projections), $(num_blocks), itmax=100_000, EPSVAL=1e-5)
            print_results(method, num_rows, num_cols, xSol, num_blocks, results, elap_time, c, df_results)
        end
    end
    timenow = Dates.now()
    filename = savename("BBGS24_Synthetic_Affine_Subspaces_parallel", (@dict timenow), "csv")
    CSV.write(datadir("sims",filename), df_results)
    
    @test true
end

##
filename = "BBGS24_Synthetic_Affine_Subspaces_parallel_timenow=2025-04-30T14:29:44.750.csv"
df_results = CSV.read(datadir("sims", filename), DataFrame)
sort!(df_results,[:num_blocks])
gdf = groupby(df_results, [:method, :num_blocks, :num_rows, :num_cols])
summary_filename = "BBGS24_Synthetic_Affine_Subspaces_parallel_summary_timenow=2025-04-30T14:29:44.750.csv"
CSV.write(datadir("sims", summary_filename), combine(gdf, :wall_clock_time => mean, :total_projections => mean))
    


##

# using MAT
# using SparseArrays
# @testset "Shepp-Logan Phantom" begin
    # filename = datadir("exp_raw/phantom.mat")
    # mf = matread(filename)
    # A = sparse(mf["A"])
    # J = Array{Int64}(undef, 0)
    # for row in eachrow(A)  xSol = A'* w
      
    #     if ~iszero(row)
    #         push!(J, row.indices[1])
    #     end
    # end
    # A = A[J, :] 
    # b = (mf["b"])[J]
#     xSol = (mf["x"])[:]
#     Xkacz = (mf["Xkacz"])[:]
#     Xsymk = (mf["Xsymk"])[:]
#     Xrand = (mf["Xrand"])[:]
#     x₀ = zeros(size(A, 2))
#     for num_blocks in [4, 16, 32]
#         @info "Running CRM with $num_blocks blocks"
#         Projections = projection_block_QR(A, b, num_blocks)
#         results = SPM(x₀, Projections, num_blocks, itmax=1000)
#         elap_time = @elapsed CRM(x₀, Projections, num_blocks, itmax=10000)
#         print_results(:CRM, xSol, num_blocks, results, elap_time)
#     end
# end




# ##

# @info "Using CUDA"
# @info CUDA.versioninfo()
# # @testset "Second Example: Random Matrix in GPU" begin
# Random.seed!(42)
# T = Float32
# num_rows, num_cols = 5500, 10_000
# A = randn(T, num_rows, num_cols)
# w = randn(T, num_rows)
# xSol = A' * w
# normalize!(xSol)
# b = A * xSol
# x₀ = zeros(T, num_cols)
# Ablock = cu(A)
# bblock = cu(b)
# xstar = cu(x₀)
# projAgpu = proj_factory_QR(Ablock, bblock)
# @benchmark CUDA.@sync projAgpu($xstar)
# projAcpu = proj_factory_QR(A, b)
# @benchmark projAcpu($x₀)

# num_blocks = 200
# ProjectionGPU = projection_block_QR(Ablock, bblock, num_blocks)
# ProjectionCPU = projection_block_QR(A, b, num_blocks)
# weights = ones(T, num_blocks) ./ num_blocks
# weightsgpu = cu(weights) 

# @benchmark CUDA.@sync mapreduce((proj) -> proj($xstar), +, $ProjectionGPU)
# @benchmark mapreduce((proj) -> proj($x₀), +, $ProjectionCPU)

##

# using SuiteSparseMatrixCollection
# using SparseArrays
# using MatrixMarket
# ssmc = ssmc_db()
# paths =  fetch_ssmc(ssmc_matrices(ssmc, "", "Franz1"), format="MM")
# paths = vcat(paths, fetch_ssmc(ssmc_matrices(ssmc, "Meszaros", "crew1"), format="MM"))
# downloaded_files = installed_ssmc()

# @testset "Matrices  from Matrix Market" begin

#     for path in paths
#         matrix_name = readdir(path)[1]
#         T = Float64
#         M = float.(mmread(path * "/" * matrix_name))
#         num_rows, num_cols = size(M)
#         if num_cols > num_rows
#             M = permutedims(M)
#             num_rows, num_cols = size(M)
#         end
#         w = randn(T, num_rows)
#         xSol = M' * w
#         normalize!(xSol)
#         b = M * xSol
#         x₀ = zeros(T, num_cols)
#         # Ablock = M[1:112,:]
#         # bblock = b[1:112]
#         # SF = qr(Ablock', ColumnNorm())
#         # projA = proj_factory_Krylov(Ablock, bblock)
#         # projB = proj_factory_QR(Ablock, bblock, SF)
#         # projC = proj_factory_QRMumps(Ablock, bblock)
#         num_blocks_vec = 20:10:100
#         @info "Solving system $(matrix_name) from matrix Market with size $(size(M)) and using several projection factories"

#         for proj in [:projection_block_QR], num_blocks in num_blocks_vec, method in [:pCRM, :CRM] 
#             @info "Running $(method) with $num_blocks blocks and $(proj) projection factory"
#             proj_bloc = eval(proj)
#             Projections = proj_bloc(M, b, num_blocks)
#             func = eval(method)
#             results = nothing
#             try
#                 results = func(x₀, Projections, num_blocks, itmax=100_000, EPSVAL=1e-5)
#             catch e
#                 @warn "Error $e in $(method) with $num_blocks blocks and $(proj) projection factory"
#                 continue
#             end
#             elap_time = @belapsed $(func)($(x₀), $(Projections), $(num_blocks), itmax=100_000)
#             !isnothing(results) && print_results(method, xSol, num_blocks, results, elap_time)
#         end
#     end
# end
# k = 1
# for proj in Projections
#     println("Projection $k")
#     proj(x₀)
#     k += 1
    
# end

# proj = Projections[6]
# proj(x₀)

# size_block = div(num_rows,  num_blocks)
