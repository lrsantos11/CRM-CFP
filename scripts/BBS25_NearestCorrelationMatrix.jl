##################################################################
### Methods
##################################################################

"""
        Boyle_Dykstra_NCM(A::Matrix{T}, ProjS, ProjU) where {T<:Real}

    Boyle-Dykstra's algorithm for finding the Nearest Correlation matrix to a given matrix A
    A is a matrix
    ProjS is the projection onto the set of positive semidefinite matrices
    ProjU is the projection onto the set of correlation matrices
    max_iter is the maximum number of iterations
    tol is the tolerance for stopping the iterations
    returns the nearest correlation matrix to A
    Reference: Higham, Nicholas J. "Computing the nearest correlation matrix—a problem from finance." IMA Journal of Numerical Analysis 22.3 (2002): 329-343.

"""
function Boyle_Dykstra_NCM(A::Matrix{T},
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


function BestCirc_V1(X₀::Matrix{T},
                  ProjS,
                  ProjU;
                  max_iter::Int64=1000, tol::Float64=1e-7) where {T<:Real}
    
    X₀ = ProjU(X₀)
    ReflectU = x -> Reflection(x, ProjU)
    Xₖ = CRMiteration(X₀, ProjS(X₀), ReflectU)
    ProjS_Xₖ = ProjS(Xₖ)

    solved = false
    tired = false
    iter = 0
    while !(solved || tired)
        Xₖ_ProjS_Xₖ = Xₖ - ProjS_Xₖ
        X₀_Xₖ = X₀ - Xₖ

        # Update Xₖ as the Projection of X₀ onto  Wₖ ∩ Hₖ ∩ U
        ProjectHₖ = x -> ProjHyperplane(Xₖ_ProjS_Xₖ, dot(Xₖ_ProjS_Xₖ, ProjS_Xₖ), x)
        ProjectWₖ = x -> ProjHyperplane(X₀_Xₖ, dot(X₀_Xₖ, Xₖ), x)
        Xₖ = Project_W_H_U(X₀, ProjU, ProjectWₖ, ProjectHₖ) 

        ProjS_Xₖ = ProjS(Xₖ)
        iter += 1
        tired = iter > max_iter
        tolBestCirc = norm(Xₖ - ProjS_Xₖ)
        solved = tolBestCirc < tol 
        @info "iter = $iter, tolBestCirc = $tolBestCirc"
    end
    return Xₖ, iter
end



"""
    BestCirc(A::Matrix{T}, ProjS, ProjU) where {T<:Real}
    Using the Best Approximation circumcentered-reflection method to find the nearest correlation matrix to A
"""


function BestCirc_Old_equiv_CRM(X₀::Matrix{T},
    ProjS,
    ProjU;
    max_iter::Int64=1000, tol::Float64=1e-8) where {T<:Real}

    X₀ = ProjU(X₀)
    Xₖ = copy(X₀)

    solved = false
    tired = false
    iter = 0
    ReflectU = x -> Reflection(x, ProjU)
    ProjS_Xₖ = ProjS(Xₖ)
    while !(solved || tired)
        Wₖ = CRMiteration(Xₖ, ProjS_Xₖ, ReflectU)
        if iter == 0
            Xₖ = Wₖ
        else
            Xₖ_ProjS_Xₖ = Xₖ - ProjS_Xₖ
            X₀_Wₖ = X₀ - Wₖ

            # Update Xₖ as the Projection of X₀ onto  Wₖ ∩ Hₖ ∩ U
            ProjectHₖ = x -> ProjHyperplane(Xₖ_ProjS_Xₖ, dot(Xₖ_ProjS_Xₖ, ProjS_Xₖ), x)
            ProjectWₖ = x -> ProjHyperplane(X₀_Wₖ, dot(X₀_Wₖ, Wₖ), x)
            Xₖ = Project_W_H_U(X₀, ProjU, ProjectWₖ, ProjectHₖ)

        end
        ProjS_Xₖ = ProjS(Xₖ)
        iter += 1
        tired = iter > max_iter
        tolBestCirc = norm(Xₖ - ProjS_Xₖ)
        solved = tolBestCirc < tol
        @info "iter = $iter, tolBestCirc = $tolBestCirc"
    end
    return Xₖ, iter
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
              ) = y +  ((b - dot(a,y))/dot(a,a)) * a



## Projection onto the intersecion of U with two halfspaces W and H
function Project_W_H_U(X₀, ProjU, ProjectW, ProjectH)

    # Projection onto the intersection of U with W
    Y = Project_U_cap_Hyperplane(X₀, ProjU, ProjectW)

    #  Check if Y is in H. If so, it is the projection
    if ProjectH(Y) ≈ Y
        return Y
    end

    Project_U_cap_H = x -> Project_U_cap_Hyperplane(x, ProjU, ProjectH)
    Reflect_U_cap_H = x -> Reflection(x, Project_U_cap_H)
    ProjU_H_X₀ = Project_U_cap_H(X₀)
    # Return the circumcenter of X₀, R_W(X₀) and R_HR_W(X₀)
    Proj_U_W_X₀ = CRMiteration(X₀, Y, Reflect_U_cap_H)
    vec_dist =  [X₀ - ProjU_H_X₀,
        X₀ - Proj_U_W_X₀]
    _, min_index = findmin(norm.(vec_dist))
    if min_index == 1
        return ProjU_H_X₀
    else
        return Proj_U_W_X₀
    end
end

function Project_U_cap_Hyperplane(X, ProjU, ProjHyp)
    ReflectH = x -> Reflection(x, ProjHyp)
    ReflectU = x -> Reflection(x, ProjU)
    return CRMiteration(X, ReflectH, ReflectU)
end


##################################################################
### Examples
##################################################################


## Example 0 from Higham (2002)
A = Float64[1 1 0
            1 1 1
            0 1 1]
# Solving with Boyle-Dykstra's algorithm

X_BD, iter = Boyle_Dykstra_NCM(A, ProjS, ProjU, max_iter=1000, tol=1e-8)

norm(X_BD - A)

q = [-0.4814, 0.7324, -0.4814]
X_BD * q

# Solving with CRM
ResutlsCRM = CRM(A, ProjS, ProjU)
X_CRM = ResutlsCRM.xApprox

norm(X_CRM - A)
norm(X_CRM - X_BD)


# Solving with BestCirc
X_BC, iter = BestCirc_V1(A, ProjS, ProjU)

norm(X_BC  - A)
X_BC * q


# Solving with BestCirc_Old
X_BC_old, iter = BestCirc_Old_equiv_CRM(A, ProjS, ProjU)

norm(X_BC_old - A)
X_BC * q




# Example 1 from Higham (2002)

A = diagm(0 => [2, 2, 2, 2.0], -1 => [-1, -1, -1.0], 1 => [-1, -1, -1.0])
X_BD, iter = Boyle_Dykstra_NCM(ProjU(A), ProjS, ProjU)

X_BD
norm(A - X_BD)

ResutlsCRM = CRM(A, ProjS, ProjU)
X_CRM = ResutlsCRM.xApprox
norm(X_CRM - A) 



X_BC, iter = BestCirc_V1(A, ProjS, ProjU)

norm(X_BC  - A)

