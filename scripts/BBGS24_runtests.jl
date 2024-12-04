using DrWatson
@quickactivate Symbol("CRM-CFP")

using Test
# Include the CRM-CFP.jl file
include(srcdir("CRM-CFP.jl"))

function print_results(method, xSol, num_blocks, results, elap_time)
    @info "Distance to solution:  $(norm(xSol - results.xApprox))"
    @info "Total Projections:  $(results.iter_total)"
    @info "Wall-clock time for $(method) with $num_blocks blocks: $(elap_time)"
    println("="^80)
end

##
@testset "First Example: Random Matrix" begin
    Random.seed!(42)
    T = Float64
    num_rows, num_cols = 10_000, 1_000
    A = randn(T, num_rows, num_cols)
    w = randn(T, num_rows)
    xSol = A'* w
    normalize!(xSol)
    b = A * xSol
    x₀ = zeros(T, num_cols)
    num_blocks_vec = 15:20
    @info "Solving system with size $(size(A)) and using Julia's QR"
    println("="^80)
    println("="^80)
    for  num_blocks in num_blocks_vec, method in [:pCRM, :CRM] 
        @info "Running $(method) with $num_blocks blocks"
        Projections = projection_block_QR(A, b, num_blocks)
        func = eval(method)
        results = func(x₀, Projections, num_blocks, itmax=100_000, EPSVAL=1e-4)
        elap_time = @belapsed $(func)($(x₀), $(Projections), $(num_blocks), itmax=100_000)
        print_results(method, xSol, num_blocks, results, elap_time)
    end
    
    @test true
end


##

# using MAT
# using SparseArrays
# @testset "Shepp-Logan Phantom" begin
    # filename = datadir("exp_raw/phantom.mat")
    # mf = matread(filename)
    # A = sparse(mf["A"])
    # J = Array{Int64}(undef, 0)
    # for row in eachrow(A)
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

using SuiteSparseMatrixCollection
using SparseArrays
using MatrixMarket
ssmc = ssmc_db()
paths = fetch_ssmc(ssmc_matrices(ssmc, "", "Franz1"), format="MM")
downloaded_files = installed_ssmc()
T = Float64
M = float.(mmread(paths[1] * "/Franz1.mtx"))
num_rows, num_cols = size(M)
w = randn(T, num_rows)
xSol = M' * w
normalize!(xSol)
b = M * xSol
x₀ = zeros(T, num_cols)
Ablock = M[1:112,:]
bblock = b[1:112]
SF = qr(Ablock', ColumnNorm())
projA = proj_factory_Krylov(Ablock, bblock)
projB = proj_factory_QR(Ablock, bblock, SF)
projC = proj_factory_QRMumps(Ablock, bblock)

@test norm(projA(x₀) - projB(x₀)) < 1e-8
@test projC(x₀) ≈ projB(x₀)

@benchmark projA($x₀)
@benchmark projB($x₀)
@benchmark projC($x₀)



num_blocks = 50
ProjectionsKrylov = projection_block_Krylov(M, b, num_blocks)
ProjectionsQR = projection_block_QR(M, b, num_blocks)
ProjectionsQRMumps = projection_block_QRMumps(M, b, num_blocks)

@benchmark ThreadsX.mapreduce((proj) -> proj($x₀), +, $ProjectionsKrylov)
@benchmark mapreduce((proj) -> proj($x₀), +, $ProjectionsQR)
@benchmark mapreduce((proj) -> proj($x₀), +, $ProjectionsQRMumps)``
@benchmark ThreadsX.mapreduce((proj) -> proj($x₀), +, $ProjectionsQRMumps)
@benchmark _test1([$x₀], $ProjectionsQRMumps)
