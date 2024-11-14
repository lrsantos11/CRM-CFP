using DrWatson
@quickactivate "CRM-CFP"

using Test
# Include the CRM-CFP.jl file
include(srcdir("CRM-CFP.jl"))


@testset "First Example: Random Matrix" begin
    Random.seed!(42)
    T = Float64
    num_rows, num_cols = 1000, 100_000
    A = randn(T, num_rows, num_cols)
    w = randn(T, num_rows)
    xSol = A'* w
    normalize!(xSol)
    b = A * xSol
    x₀ = zeros(T, num_cols)
    for method in [:SPM, :CSRM], num_blocks in [5, 10, 20]
        @info "Running $(method) with $num_blocks blocks"
        Projections = simultaneousproj_IndAffine(A, b, num_blocks)
        func = eval(method)
        results = func(x₀, Projections, num_blocks, itmax=10000)
        @info norm(xSol - results.xApprox)
        @info results.iter_total
        @time func(x₀, Projections, num_blocks, itmax=10000)
        # @info "Wall-clock time for SPM with $num_blocks blocks: $t"
    end
    for num_blocks in [5, 10, 20]
        @info "Running CRM with $num_blocks blocks"
        Projections = successiveproj_IndAffine(A, b, num_blocks)
        results = CRM(x₀, Projections, num_blocks, itmax=10000)
        @info norm(xSol - results.xApprox)
        @info results.iter_total
        @time CRM(x₀, Projections, num_blocks, itmax=10000)
        # @info "Wall-clock time for SPM with $num_blocks blocks: $t"
    end
    
    
    @test true
end

