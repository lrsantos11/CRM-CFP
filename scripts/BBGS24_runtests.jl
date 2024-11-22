using DrWatson
@quickactivate "CRM-CFP"

using Test
# Include the CRM-CFP.jl file
include(srcdir("CRM-CFP.jl"))

function print_results(method, xSol, num_blocks, results, elap_time)
    @info "Distance to solution:  $(norm(xSol - results.xApprox))"
    @info "Total Projections:  $(results.iter_total)"
    @info "Wall-clock time for $(method) with $num_blocks blocks: $(elap_time)"
    println("="^80)
end

@testset "First Example: Random Matrix" begin
    Random.seed!(42)
    T = Float64
    num_rows, num_cols = 500, 10_000
    A = randn(T, num_rows, num_cols)
    w = randn(T, num_rows)
    xSol = A'* w
    normalize!(xSol)
    b = A * xSol
    x₀ = zeros(T, num_cols)
    block_size = [20, 50, 100]
    @info "Solving system with size $(size(A)) and using Julia's QR"
    println("="^80)
    println("="^80)
    for method in [:CSRM, :SPM, :CRM],  num_blocks in block_size
        @info "Running $(method) with $num_blocks blocks"
        Projections = projection_block_QR(A, b, num_blocks)
        func = eval(method)
        results = func(x₀, Projections, num_blocks, itmax=50000, EPSVAL=1e-4)
        elap_time = @belapsed $(func)($(x₀), $(Projections), $(num_blocks), itmax=10000)
        print_results(method, xSol, num_blocks, results, elap_time)
    end
    
    @test true
end







##

# using MAT
# using SparseArrays
# @testset "Shepp-Logan Phantom" begin
#     filename = datadir("exp_raw/phantom.mat")
#     mf = matread(filename)
#     A = sparse(mf["A"])
#     J = Array{Int64}(undef, 0)
#     for row in eachrow(A)
#         if ~iszero(row)
#             push!(J, row.indices[1])
#         end
#     end
#     A = A[J, :] 
#     b = (mf["b"])[J]
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


