##################################################################
## Basic Functions for Ellipsoids tests and plots
##################################################################
using CSV
using DataFrames
using BenchmarkProfiles
using BenchmarkTools
using LazySets
import LazySets: Ellipsoid
import Base: in
"""
Structure of an Ellipsoid satisfying dot(x,A*x) + 2*dot(b,x) ≤ α
"""
@kwdef struct EllipsoidCRM
    A::AbstractMatrix
    b::AbstractVector
    α::Number
end

function in(x₀::Vector, ell::EllipsoidCRM)
    A, b, α = ell.A, ell.b, ell.α
    if dot(x₀, A * x₀) + 2 * dot(b, x₀) ≤ α
        return true
    else
        return false
    end
end
"""
Transform Ellipsoid in format dot(x-c,Q⁻¹*(x-c)) ≤ 1 
into format  dot(x,A*x) + 2*dot(b,x) ≤ α
from shape matrix Q and center of ellipsoid c
"""
function EllipsoidCRM(c::Vector, Q::AbstractMatrix)
    A = inv(Matrix(Q))
    b = -A * c
    α = 1 + dot(c, b)
    return EllipsoidCRM(A, b, α)
end


"""
EllipsoidCRM(ell)

Transform Ellipsoid in format dot(x-c,Q⁻¹*(x-c)) ≤ 1 from LazySets
into format  dot(x,A*x) + 2*dot(b,x) ≤ α
from shape matrix Q and center of ellipsoid c
"""
EllipsoidCRM(ell::Ellipsoid) = EllipsoidCRM(ell.center, ell.shape_matrix)


"""
Ellipsoid(ell)
Transform Ellipsoid in format  dot(x,A*x) + 2*dot(b,x) ≤ α to format dot(x-c,Q⁻¹*(x-c)) ≤ 1 from LazySets
from shape matrix Q and center of ellipsoid c
"""
function Ellipsoid(ell::EllipsoidCRM)
    c = -(ell.A \ ell.b)
    β = ell.α - dot(c, ell.b)
    Q = Symmetric(inv(Matrix(ell.A) / β))
    return Ellipsoid(c, Q)
end



"""
Proj_Ellipsoid(x₀, ell)
Projects x₀ onto ell, and EllipsoidCRM using an ADMM algorithm as reported by Jia, Cai and Han [Jia2007]

[Jia2007] Z. Jia, X. Cai, e D. Han, “Comparison of several fast algorithms for projection onto an ellipsoid”, Journal of Computational and Applied Mathematics, vol. 319, p. 320–337, ago. 2017, doi: 10.1016/j.cam.2017.01.008.
"""
function Proj_Ellipsoid(x₀::Vector,
    ell::EllipsoidCRM;
    itmax::Int=10_000,
    ε::Real=1e-8,
    verbose::Bool=false)
    x₀ ∉ ell ? nothing : return x₀
    A, b, α = ell.A, ell.b, ell.α
    ϑₖ = 10 / norm(A)
    n = length(x₀)
    B = sqrt(Matrix(A))
    issymmetric(B) ? BT = B : BT = B'
    b̄ = B \ (-b)
    αplusb̄2 = α + norm(b̄)^2
    r = sqrt(αplusb̄2)
    yₖ = ones(n)
    λₖ = ones(n)
    xₖ = x₀
    it = 0
    tolADMM = 1.0
    function ProjY(y)
        normy = norm(y)
        if αplusb̄2 - normy^2 ≥ 0.0
            return y
        else
            return (r / normy) * y
        end
    end
    Ā = (I + ϑₖ * A)
    normRxₖ = 0.0
    while tolADMM ≥ ε^2 && it ≤ itmax
        uₖ = x₀ + BT * (λₖ + ϑₖ * (yₖ + b̄))
        xₖ = Ā \ uₖ
        wₖ = B * xₖ - λₖ ./ ϑₖ - b̄
        normwₖ = norm(wₖ)
        normwₖ ≤ r ? yₖ = wₖ : yₖ = (r / normwₖ) * wₖ
        Rxₖ = xₖ - x₀ - BT * λₖ
        Ryₖ = yₖ - ProjY(yₖ - λₖ)
        Rλₖ = B * xₖ - yₖ - b̄
        λₖ -= ϑₖ * Rλₖ
        normRxₖ = norm(Rxₖ)
        normRyₖ = norm(Ryₖ)
        normRλₖ = norm(Rλₖ)
        tolADMM = sum([normRxₖ^2, normRyₖ^2, normRλₖ^2])
        it += 1
        # if normRxₖ < normRλₖ*(0.1/n)
        #     ϑₖ *= 2
        # elseif normRxₖ > normRλₖ*(0.9/n)
        #     ϑₖ *= 0.5
        # end

    end

    verbose && @info it, normRxₖ
    return real.(xₖ)
end



"""
Function value of ellipsoid
"""
function func_EllipsoidCRM(x::Vector, ell::EllipsoidCRM)
    A, b, α = ell.A, ell.b, ell.α
    return dot(x, A * x) + 2 * dot(b, x) - α
end


"""
Ellipsoid initial point 
"""
function InitalPoint_EllipsoidCRM(Ellipsoids::Vector{EllipsoidCRM}, n::Int; ρ::Number=1.2)
    x₀ = StartingPoint(n)
    iter_starting_point = 1
    while any(Ref(x₀) .∈ Ellipsoids) && iter_starting_point < 100
        iter_starting_point += 1
        x₀ .*= ρ
    end
    return x₀
end

"""
    Approximate projection onto ellipsoid using the (perturbed) Separating Hyperplane defined using subdiferential ``u``of function ``f``.
    ```math
    P_{S_i^k}(x^k) = x^k -\\frac{\\max\\{0,f_i(x^k)+1\\epsilon_k\\}} {\\|u \\|^2}u,
    ```
"""


function ApproxProj_Ellipsoid(x::Vector,
                              Ellipsoid::EllipsoidCRM;
                              λ::Real = 1.0, # relaxation parameter
                              ϵ::Real = 0.0 # perturbation
                              )
    A, b = Ellipsoid.A, Ellipsoid.b
    Ax = A * x
    fx = func_EllipsoidCRM(x, Ellipsoid)
    ∂fx = 2 * (Ax + b)
    
    return λ * (x .- (max(0.0,fx + ϵ) / dot(∂fx, ∂fx)) * ∂fx) .+ (1 - λ) * x
end
    
    """
        Approximate projection onto ellipsoid given as a dict 
    """

function  ApproxProj_Ellipsoid(x::Vector,Ellipsoid::Dict; kwargs...)
    @unpack A, b, α  = Ellipsoid
    ell = EllipsoidCRM(A,b,α)
    return ApproxProj_Ellipsoid(x,ell; kwargs...)
end





"""
    Approximate projection onto ellipsoid using ProdSpace
"""

function ApproxProjectEllipsoids_ProdSpace(X::Vector,
    Ellipsoids::Vector{EllipsoidCRM})
    proj = similar(X)
    for index in eachindex(proj)
        proj[index] = ApproxProj_Ellipsoid(X[index], Ellipsoids[index])
    end
    return proj
end

##################################################################
### Creating Ellipsoids Samples
##################################################################

"""
SampleTwoEllipsoids(n, p; λ=1.1, γ=1.5)
Creates two ellipsoids in ℝⁿ that intersect. The intersection is regulated by λ.

"""
function SampleTwoEllipsoids(n::Int,  # dimension
    p::Real; # sparsity of matrix A
    λ::Real=1.1, # parameter for touching ellipsoid
    γ::Real=1.5 # parameter for making A positive definite
)
    Ellipsoids = EllipsoidCRM[]
    A = Matrix(sprandn(n, n, p))
    A = (γ * I + A' * A)
    a = rand(n)
    b = A * a
    adotAa = dot(a, b)
    b .*= -1.0
    α = (1 + γ) * adotAa
    push!(Ellipsoids, EllipsoidCRM(A, b, α))
    TouchEll, Center2, TouchPoint = GenerateTouchingEllipsoid(Ellipsoids[1], n, λ=λ)
    push!(Ellipsoids, TouchEll)
    return Ellipsoids, Center2, TouchPoint
end




"""
SampleEllipsoids(n, m, p; λ=1.1, γ=1.5)
Creates m ellipsoids in ℝⁿ that intersect. The intersection is regulated by λ.

"""
function SampleEllipsoids(n::Int,  # dimension
    m::Int,  # number of ellipsoids
    p::Real; # sparsity of matrix A
    λ::Real=1.1, # parameter for touching ellipsoid
    γ::Real=1.5 # parameter for making A positive definite
)
    Ellipsoids, CenterEll2, TouchPoint = SampleTwoEllipsoids(n, p, λ=λ, γ=γ)
    point_inter = 0.5 * ((1 - λ)CenterEll2 + (1 + λ)TouchPoint)
    for _ in 3:m
        center = randn(n)
        while any(Ref(center) .∈ Ellipsoids)
            center *= 1.5
        end
        d = λ * (point_inter - center)
        push!(Ellipsoids, GenerateEllipsoid(center, d))
    end
    return Ellipsoids, point_inter
end

"""
    GenerateTouchingEllipsoid(ell, n; λ=1.0)
    Given an ellipsoid ell, generates a touching ellipsoid in ℝⁿ

    """


function GenerateTouchingEllipsoid(ell::EllipsoidCRM,
    n::Int;
    λ::Real=1.1)
    c = randn(n)
    while c ∈ ell
        c *= 1.5
    end
    x̂ = Proj_Ellipsoid(c, ell)
    d = λ * (x̂ - c)
    return GenerateEllipsoid(c, d), c, x̂
end



"""
 GenerateEllipsoid(center, semi_axis)
 Given center and larger semi_axis, generates an ellipsoid in ℝⁿ
"""

function GenerateEllipsoid(center::Vector,
    semi_axis::Vector)
    n = length(center)
    semi_axis_norm = norm(semi_axis)
    Λ = Diagonal([semi_axis_norm; semi_axis_norm .+ 2 * rand(n - 1)])
    # Λ = Diagonal([semi_axis_norm; rand(n - 1)* .8 * semi_axis_norm])   
    Q, _ = qr(randn(n, n))
    M2 = Q' * Λ .^ 2 * Q
    return EllipsoidCRM(center, 0.5(M2 + M2'))
end


