##################################################################
### Methods
##################################################################

"""
        Dykstra(A::Matrix{T}, ProjS, ProjU) where {T<:Real}

    Boyle-Dykstra's algorithm for finding the Nearest Correlation matrix to a given matrix A
    A is a matrix
    ProjS is the projection onto the set of positive semidefinite matrices
    ProjU is the projection onto the set of correlation matrices
    max_iter is the maximum number of iterations
    tol is the tolerance for stopping the iterations
    returns the nearest correlation matrix to A
    Reference: Higham, Nicholas J. "Computing the nearest correlation matrix—a problem from finance." IMA Journal of Numerical Analysis 22.3 (2002): 329-343.

"""
function Boyle_Dykstra(A::Matrix{T},
                        ProjS,
                        ProjU;
                        max_iter::Int64=1000, tol::Float64=1e-8) where {T<:Real}
    
    num_rows, num_cols = size(A)
    ΔS = zeros(T, num_rows, num_cols)
    Y = copy(A)
    solved = false
    tired = false
    iter = 0
    while !(solved || tired)
        R = Y - ΔS
        X = ProjS(R)
        ΔS = X - R
        Yold = copy(Y)
        Y .= ProjU(X)
        iter += 1
        iter > max_iter && (tired = true)
        norm(Y - Yold) < tol && (solved = true)
    end
    return Y, iter
end


"""
    BestCirc(A::Matrix{T}, ProjS, ProjU) where {T<:Real}
    Using the Best Approximation circumcentered-reflection method to find the nearest correlation matrix to A
"""


function BestCirc(X₀::Matrix{T},
                  ProjS,
                  ProjU;
                  max_iter::Int64=1000, tol::Float64=1e-8) where {T<:Real}
    
    X₀ = ProjU(X₀)
    Xₖ = copy(X₀)

    solved = false
    tired = false
    iter = 0
    ReflectS(x) = Reflection(x, ProjS)
    ReflectU(x) = Reflection(x, ProjU)
    ProjS_Xₖ = ProjS(Xₖ)
    while !(solved || tired)
        XₖOld = copy(Xₖ)
        Xₖ_PSXₖ = Xₖ - ProjS_Xₖ
        ProjectHₖ(x) = ProjHyperplane(Xₖ_PSXₖ, dot(Xₖ_PSXₖ, ProjS_Xₖ), x)
        ReflectHₖ(x) = Reflection(x, ProjectHₖ)
        Wₖ = CRMiteration(Xₖ, ProjS_Xₖ, ReflectU)
        X₀_Wₖ = X₀ - Wₖ
        ProjectWₖ(x) = ProjHyperplane(X₀_Wₖ, dot(X₀_Wₖ, Wₖ), x)
        ReflectWₖ(x) = Reflection(x, ProjectWₖ)

        
        Yₖ = ReflectHₖ(X₀)
        Zₖ = ReflectWₖ(Yₖ)
        Uₖ = ReflectU(Zₖ)
    
        Xₖ = FindCircumcentermSet([X₀, Yₖ, Zₖ, Uₖ])
        ProjU(Xₖ)
        iter += 1
        iter > max_iter && (tired = true)
        tolBestCirc = norm(Xₖ - XₖOld)
        tolBestCirc < tol && (solved = true)
    end
    return Results(iter_total=iter, final_tol=tolBestCirc, xApprox=Xₖ, method=BestCirc)
end




## Projection onto Semidefinite cone using ProximalOperators.jl
function ProjS(X; SemiDefiniteCone=IndPSD())
    proj, _ = prox(SemiDefiniteCone, X)
    return proj
end

## Projection onto Correlation affine subspace
function ProjU(X)
    T = eltype(X)
    proj = copy(X)
    proj[diagind(X)] .= one(T)
    return proj
end





ProjHyperplane(a::AbstractArray, # a is the orthogonal vector
              b::Number, # right-hand side value 
              y::AbstractArray # point to be projected
              ) = y +  ((b - dot(a,y))/norm(a,2)) * a




##################################################################
### Examples
##################################################################


# Example 1 from Higham (2002)
A = Float64[1 1 0
            1 1 1
            0 1 1]
# Solving with Boyle-Dykstra's algorithm
X, iter = Boyle_Dykstra(A, ProjS, ProjU)

norm(corr_A - A)

q = [-0.4814, 0.7324, -0.4814]
X*q


"""
    CRM(A::Matrix{T}, ProjS, ProjU) where {T<:Real}
    Using the Circumcentered-Reflection method to find the nearest correlation matrix to A
"""
resultsCRM = CRM(A, ProjU, ProjS)

xCRM = resultsCRM.xApprox

norm(X - xCRM)