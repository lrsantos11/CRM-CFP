__precompile__()

using DelimitedFiles
using Dates
using Distributions
using LinearAlgebra
using GenericLinearAlgebra
using PolynomialRoots
using Printf
using ProximalOperators
using Random
using SparseArrays

import Base.@kwdef
####################################
@kwdef mutable struct Results
    iter_total::Int
    proj_total::Int = 0
    final_tol::Number
    xApprox::AbstractVector
    method::Symbol
    date::DateTime = Dates.now()
end

Results() = Results(iter_total=0, final_tol=0.0, xApprox=[], method = :None)
####################################
"""
find_circumcenter!(CC, X)

Finds the Circumcenter of linearly independent vectors ``x_0,x_1,…,x_m``, belonging to ``X``,
as described in [^Behling2018a] and [^Behling2018b].

[^Behling2018a]: Behling, R., Bello Cruz, J.Y., Santos, L.-R.:
Circumcentering the Douglas–Rachford method. Numer. Algorithms. 78(3), 759–776 (2018).
[doi:10.1007/s11075-017-0399-5](https://doi.org/10.1007/s11075-017-0399-5)
[^Behling2018b]: Behling, R., Bello Cruz, J.Y., Santos, L.-R.:
On the linear convergence of the circumcentered-reflection method. Oper. Res. Lett. 46(2), 159-162 (2018).
[doi:10.1016/j.orl.2017.11.018](https://doi.org/10.1016/j.orl.2017.11.018)
"""
function find_circumcenter!(CC, X)
    # Finds the Circumcenter of  linearly independent points  X = [X1, X2, X3, ... Xn]
    # println(typeof(X))
    T = eltype(CC)
    lengthX = length(X)
    if lengthX == 1
        return X[1]
    elseif lengthX == 2
        return 0.5 * (X[1] + X[2])
    end
    V = []
    b = T[]
    # Forms V = [X[2] - X[1] ... X[n]-X[1]]
    # and b = [dot(V[1],V[1]) ... dot(V[n-1],V[n-1])]
    for ind in 2:lengthX
        difXnX1 = X[ind] - X[1]
        push!(V, difXnX1)
        push!(b, dot(difXnX1, difXnX1))
    end

    # Forms Gram Matrix
    dimG = lengthX - 1
    G = diagm(b)

    for irow in 1:(dimG-1)
        for icol in (irow+1):dimG
            G[irow, icol] = dot(V[irow], V[icol])
            G[icol, irow] = G[irow, icol]
        end
    end
    if isposdef(G)
        L = cholesky(G)
        y = L \ b
    else
        @warn "Circumcenter matrix is not positive definite. Circumcenter is not unique"
        y = G \ b
    end
    copyto!(CC, X[1])
    for ind in 1:dimG
        CC += 0.5 * y[ind] * V[ind]
    end
    return CC
end

FindCircumcentermSet(X) = find_circumcenter!(similar(X[1]), X)

####################################
function find_circumcenter(X::Vector{Vector{BigFloat}})
    lengthX = length(X)
    lengthX > 3 && error("Only computes BigFloat Circumcenter up to three")
    if lengthX == 1
        return X[1]
    elseif lengthX == 2
        return 0.5 * (X[1] + X[2])
    end
    x, y, z = X
    Sᵤ = y - x
    Sᵥ = z - x
    norm_Sᵤ = dot(Sᵤ, Sᵤ)
    norm_Sᵥ = dot(Sᵥ, Sᵥ)
    prod = dot(Sᵤ, Sᵥ)
    A_inv = [norm_Sᵥ -prod; -prod norm_Sᵤ] ./ (norm_Sᵤ * norm_Sᵥ - prod^2)
    b = [1 / 2 .* norm_Sᵤ; 1 / 2 .* norm_Sᵥ]
    sol = A_inv * b
    C = x + sol[1] * Sᵤ + sol[2] * Sᵥ
    return C
end
####################################
"""
    proj = ProjectIndicator(indicator,x)
    Projection using Indicator Function from `ProximalOperators.jl`
    """
function ProjectIndicator(indicator, x)
    proj, fproj = prox(indicator, x)
    return proj
end

####################################
"""
reflec = ReflectIndicator(indicator,x)
Reflection using Indicator Function from `ProximalOperators.jl`
"""
function ReflectIndicator(indicator, x)
    proj = ProjectIndicator(indicator, x)
    reflec = 2 * proj - x
end
####################################
function StartingPoint(n::Int64; max_value::Number=15, min_value::Number=5)
    ## Creates a random point in R^n
    x = zeros(n)
    while norm(x) < 2
        x = randn(n)
    end
    # norm between min_value and max_value
    foonorm = (max_value - min_value) * rand() + min_value
    return foonorm * x / norm(x)

end

####################################
function StartingPoint(m::Int64, n::Int64)
    ## Creates a random point in R^{m\times n}
    X = zeros(m, n)
    for j in 1:n
        X[:, j] = StartingPoint(m)
    end
    return X
end
####################################
function GenerateTwoAffines(n::Int64, affine::Bool=true)
    ma = rand(1:n-1) ## number of extra normals of A
    mb = rand(1:n-1-ma) ## number of extra normals of B
    A = randn(ma, n)
    B = randn(mb, n)
    if affine
        a = randn(ma)
        b = randn(mb)
    else
        a = zeros(ma)
        b = zeros(mb)
    end

    return A, a, ma, B, b, mb
end
####################################

"""
AffineRandomPair(n,cmax,affine)

Generates a pair of random affine spaces with intersection 
of dimension `cmax`.
"""
function GenerateTwoIntersectingAffines(n::Int,
    cmax::Int;
    affine::Bool=false)
    cmax >= 2 || error("cmax must be greater than 3")
    n >= cmax || error("n must be greater than intersection dimension")
    @show mcex = rand(1:cmax) ## number of common normals
    maex = rand(1:n-3*mcex) ## number of extra normals of A
    mbex = rand(1:n-2*mcex-maex) ## number of extra normals of B

    Cex = randn(mcex, n)
    while !(rank(Cex) == mcex)
        Cex = randn(mcex, n)
    end

    Aex = randn(maex, n)
    while !(rank([Aex; Cex]) == maex + mcex)
        Aex = randn(maex, n)
    end

    Bex = randn(mbex, n)
    while !(rank([Aex; Cex; Bex]) == maex + mcex + mbex)
        Bex = randn(mbex, n)
    end

    ma = maex + mcex
    mb = mbex + mcex

    # Space A
    A = [Aex; Cex]

    # Space B
    B = [Cex; Bex]

    if affine
        a = randn(ma)
        b = randn(mb)
    else
        a = zeros(ma)
        b = zeros(mb)
    end

    AffineA = IndAffine(A, a)
    AffineB = IndAffine(B, b)
    Affines = [AffineA, AffineB]

    return Affines
end
####################################

function printOnFile(filename::String, printline::AbstractArray; deletefile::Bool=false)
    if isempty(filename)
        return
    end
    if deletefile
        try
            rm(filename)
        catch e
            # @warn e
        end
    end
    open(filename, "a") do file
        writedlm(file, printline)
    end
end

####################################
function printOnFile(filename::String, k::Int, tolCRM::Number, xCRM::AbstractVector; deletefile::Bool=false, isprod::Bool=false)
    if isempty(filename)
        return
    end
    if deletefile
        try
            rm(filename)
        catch e
            # @warn e
        end
    end
    isprod ? printline = hcat(k, tolCRM, xCRM[1]') : printline = hcat(k, tolCRM, xCRM')
    open(filename, "a") do file
        writedlm(file, printline)
    end
end
####################################
function FriedrichsAngleAB(A, B)
    # Calculating the Friederich Angle
    QA, RA = qr(A')
    QB, RB = qr(B')
    # S =  svd(QA'*QB)
    # Angle in Degrees
    angleAB = eigmax(QA' * QB)
    return acos(angleAB)
    # Angle in Radians
    # angleAB = acos(maximum(S))
    # println(S[2])
    # ind = findfirst(x -> x<(1-1e-8),S[2])
    # return maximum(S[2])
    # return acos(S[2][ind])
end
####################################
####################################
"""
Reflection(x,fproj)
Reflection using any projection function `fproj`
"""
function Reflection(x, fproj)
    return 2 * fproj(x) - x
end
####################################
"""
        ProjectEpiQuadratic(s,v,α)
Project ``(s,v) ∈ R^{n+1}`` onto the epigraph of ``f(x) = αx^Tx``, such that ``f(v) '≤ s``.
"""
function ProjectEpiQuadratic(s::Number, v::Union{AbstractArray,Number}, α::Number=1.0)
    if α * dot(v, v) <= s
        return s, v
    end
    #PolynomialCoefficients
    # a3μ³ + a2μ² + a1μ  + a0 = 0    
    a0 = s - α * dot(v, v)
    a1 = 4 * α * s + 1.0
    a2 = 4 * α^2 * s + 4 * α
    a3 = 4 * α^2
    r = roots([a0, a1, a2, a3])
    indexreal = findall(x -> abs.(x) < 1e-12, imag.(r))
    μ = (maximum(real.(r[indexreal])))
    x = 1 / (1 + 2 * α * μ) * v
    t = μ + s
    return t, x
end
####################################
"""
        ProjectEpiQuadratic(x,α)
Project ``x = [x₀,t₀] ∈ R^{n+1}`` onto the epigraph of ``f(u) = αu^Tu``, such that ``f(x₀) '≤ t₀   ``.
"""
function ProjectEpiQuadratic(x::AbstractArray, α::Number=1.0)
    t, x = ProjectEpiQuadratic(x[end], x[1:end-1], α)
    return [x; t]
end

####################################
"""
    ProjectBall(x, v, r)
Returns the orthogonal projection of `x` onto the L2-Ball  centered in `v` with radius `r`.
Uses the `ProximalOperators.jl` toolbox

"""
function ProjectBall(x::AbstractVector, v::AbstractVector, r::Number)
    Ball = IndBallL2(r)
    proj, fproj = prox(Translate(Ball, -v), x)
    return proj
end
####################################
"""
ProjectProdDiagonal(X,num_sets)


"""
function ProjectProdDiagonal(X::AbstractVector)
    inner_proj = mean(X)
    proj = similar(X)
    for index in eachindex(proj)
        proj[index] = inner_proj
    end
    return proj
end


####################################
"""
ProjectProdSets(X,Projections)


"""
function ProjectProdSets(X::AbstractVector, SetsProjections::AbstractVector{Function})
    proj = similar(X)
    for index in eachindex(proj)
        proj[index] = SetsProjections[index](X[index])
    end
    return proj
end
####################################

"""
ProjectSetsIndicators(X,SetsIndicators)
Projection on the Product Space of Half Spaces
"""
function ProjectSetsIndicators(X::AbstractVector, SetsIndicators::AbstractVector{T}) where {T<:Union{IndAffine,IndSOC}}
    proj = similar(X)
    for index in eachindex(SetsIndicators)
        proj[index] = ProjectIndicator(SetsIndicators[index], X[index])
    end
    return proj
end

####################################

ProjectProdSpace(X::AbstractVector, Projections::AbstractVector{Function}) = ProjectProdSets(X, Projections)
ProjectProdSpace(X::AbstractVector, Projections::AbstractVector{T}) where {T<:Union{IndAffine,IndSOC}} = ProjectSetsIndicators(X, Projections)

####################################

"""
        Tolerance(x,xold;xsol,normytpe)
"""
function Tolerance(x::AbstractVector, xold::AbstractVector, xsol::AbstractVector;
    norm_p::Number=2)
    if isempty(xsol)
        return norm(x - xold, norm_p)
    else
        return norm(x - xsol, norm_p)
    end
end

####################################
"""
    ApproxProject(x,g,∂g)

"""
function ApproxProject(x::AbstractVector, g::Function, ∂g::Function; λ::Number=1.0)
    gx = g(x)
    if gx ≤ 0
        return x
    else
        ∂gx = ∂g(x)
        return λ * (x .- (gx / dot(∂gx, ∂gx)) * ∂gx) .+ (1 - λ) * x
    end
end

"""
    createDataframes(Methods::Vector{Symbol})
    This function creates two dataframes, one for storing the results of the methods, and another for storing the filenames of the problems.
"""
function createDataframes(Methods::Vector{Symbol}; projections::Bool=false)
    dfResults = DataFrame(Problem=String[])
    dfFilenames = copy(dfResults)
    for mtd in Methods
        insertcols!(dfResults, join([mtd, "_it"]) => Int[])
        projections &&  insertcols!(dfResults, join([mtd, "_projs"]) => Int[])
        insertcols!(dfResults, join([mtd, "_tol"]) => Real[])
        insertcols!(dfResults, join([mtd, "_elapsed"]) => Real[])
        insertcols!(dfFilenames, join([mtd, "filename"]) => String[])
    end
    return dfResults, dfFilenames
end


"""
   Simultaneous Projections onto Affine Spaces
"""


function simultaneousproj_IndAffine(A::Matrix, b::Vector, num_blocks::Int)
    Indicators = []
    block_size = div(size(A, 1), num_blocks)
    for i = 1:num_blocks-1
        Ablock = A[(i-1)*block_size+1:i*block_size, :]
        bblock = b[(i-1)*block_size+1:i*block_size]
        IndAblock = IndAffine(Ablock, bblock)
        push!(Indicators, IndAblock)
    end
    Ablock = A[(num_blocks-1)*block_size+1:end, :]
    bblock = b[(num_blocks-1)*block_size+1:end]
    push!(Indicators, IndAffine(Ablock, bblock))
    function proj(x)
        map(Ind -> prox(Ind, x)[1], Indicators)
        return [prox(Ind, x)[1] for Ind in Indicators]
    end
    return proj
end


function successiveproj_IndAffine(A::AbstractMatrix, b::AbstractVector, num_blocks::Int)
    Indicators = []
    block_size = div(size(A, 1), num_blocks)
    for i = 1:num_blocks-1
        Ablock = A[(i-1)*block_size+1:i*block_size, :]
        bblock = b[(i-1)*block_size+1:i*block_size]
        IndAblock = IndAffine(Ablock, bblock, iterative = true)
        push!(Indicators, IndAblock)
    end
    Ablock = A[(num_blocks-1)*block_size+1:end, :]
    bblock = b[(num_blocks-1)*block_size+1:end]
    push!(Indicators, IndAffine(Ablock, bblock, iterative=true))
    return [x -> prox(Ind, x)[1] for Ind in Indicators]
end


####################################
using QRMumps

"""
    projection_QR(A, b, SF, xstar)
"""

function proj_factory_QR(Ablock::AbstractMatrix,
                         bblock::AbstractVector;
                         verbose::Bool=false)
    # Sparse Projection using QRMumps
    # A^T(AA^T)^{-1} = QR^{-T}
    T = eltype(Ablock)
    _, num_cols = size(Ablock)
    z_num_rows = similar(bblock)
    z_num_cols = zeros(T, num_cols)
    # Initialize QRMumps
    qrm_init()

    # Initialize data structures
    spmat = qrm_spmat_init(Ablock)
    spfct = qrm_spfct_init(spmat)

    # Specify that we want the Q-less QR factorization
    qrm_set(spfct, "qrm_keeph", 0)

    # Perform symbolic analysis of Aᵀ and factorize Aᵀ = QR
    qrm_analyse!(spmat, spfct; transp='t')
    qrm_factorize!(spmat, spfct, transp='t')

    function proj(xstar)
        verbose && @info "mul!"
        mul!(z_num_rows, Ablock, xstar) # z_num_rows = A*xstar
        verbose && @info "axpy!"
        !iszero(bblock) && axpy!(-one(T), bblock, z_num_rows) # z_num_rows = A*xstar - b
        verbose && @info "QR"
        qrm_solve!(spfct, z_num_rows, z_num_cols, transp='t')
        qrm_solve!(spfct, z_num_cols, z_num_rows, transp='n')
        mul!(z_num_cols, Ablock', z_num_rows)
        return (xstar - z_num_cols)
    end
    return proj
end

####################################

using Krylov
using LinearOperators
function proj_factory_Krylov(Ablock::AbstractMatrix,
                       bblock::AbstractVector)
    # Sparse Projection proj_factory_Krylov
    issparse(Ablock) || error("Matrix A must be sparse")
    T = eltype(Ablock)
    opA = LinearOperator(Ablock)
    cgne_solver = CgneSolver(opA, bblock)
    rhs = similar(bblock)
   
   
   function projA(xstar)
        @info "mul!"
        mul!(rhs, opA, xstar) # rhs = A*xstar
        @info "axpy!"
        !iszero(b) && axpy!(-one(T), bblock, rhs) # rhs = A*xstar - b
        @info "cgne!"
        cgne!(cgne_solver, opA, rhs, rtol=1e-10) # cgne_solver.x = Aᵀ(AAᵀ)⁻¹(A*xstar - b)
        return (xstar - cgne_solver.x)
    end
    return projA
end


