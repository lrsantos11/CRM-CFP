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
    Approximate projection onto ellipsoid
"""


function ApproxProj_Ellipsoid(x::Vector,
    Ellipsoid::EllipsoidCRM;
    λ::Real=1.0)
    A, b, α = Ellipsoid.A, Ellipsoid.b, Ellipsoid.α
    Ax = A * x
    gx = dot(x, Ax) + 2 * dot(b, x) - α
    if gx ≤ 0
        return x
    else
        ∂gx = 2 * (Ax + b)
        return λ * (x .- (gx / dot(∂gx, ∂gx)) * ∂gx) .+ (1 - λ) * x
    end
end
    
    """
        Approximate projection onto ellipsoid
    """

function  ApproxProj_Ellipsoid(x::Vector,Ellipsoid::Dict; λ::Real = 1.0)
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



