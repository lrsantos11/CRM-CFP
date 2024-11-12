##################################################################
### Methods
##################################################################



"""
        dykstra_NCM(A::Matrix{T}, ProjS, ProjU) where {T<:Real}

    Boyle-Dykstra's algorithm for finding the Nearest Correlation matrix to a given matrix A
    A is a matrix
    ProjS is the projection onto the set of positive semidefinite matrices
    ProjU is the projection onto the set of correlation matrices
    max_iter is the maximum number of iterations
    tol is the tolerance for stopping the iterations
    returns the nearest correlation matrix to A
    Reference: Higham, Nicholas J. "Computing the nearest correlation matrix—a problem from finance." IMA Journal of Numerical Analysis 22.3 (2002): 329-343.

"""
function dykstra_NCM(A::AbstractMatrix{T},
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
        Y .= ProjU(X)
        iter += 1
        iter > max_iter && (tired = true)
        toldykstra = rel_error(Y, X)
        @info "iter = $iter, toldykstra = $toldykstra"
        toldykstra < tol && (solved = true)
    end
    return Y, iter
end


"""
    bestcirc_MAP(A::Matrix{T}, ProjS, ProjU) where {T<:Real}
    Using the Best Approximation circumcentered-reflection method to find the nearest correlation matrix to A
"""


function bestcirc_MAP(X₀,
                  ProjS,
                  ProjU;
                  max_iter::Int64=1000, tol::Float64=1e-7) where {T<:Real}
    
    X̂₀ = ProjU(X₀)
    ReflectU = x -> Reflection(x, ProjU)
    Xₖ = CRMiteration(X̂₀, ProjS(X̂₀), ReflectU)
    ProjS_Xₖ = ProjS(Xₖ)


    tired = false
    iter = 1
    tolbestcirc = norm(Xₖ - ProjS_Xₖ)
    solved = tolbestcirc < tol
    @info "iter = $iter, tolbestcirc = $tolbestcirc"
    while !(solved || tired)
        Xₖ_ProjS_Xₖ = Xₖ - ProjS_Xₖ
        # X̂₀_Xₖ = X̂₀ - Xₖ
        IndHₖ = IndHalfspace(Xₖ_ProjS_Xₖ, dot(Xₖ_ProjS_Xₖ, ProjS_Xₖ))
        # IndWₖ = IndHalfspace(X̂₀_Xₖ, dot(X̂₀_Xₖ, Xₖ))
        # Note that  Xₖ is the Projection of X̂₀ onto  Wₖ ∩ U
        Xₖ = Project_W_H_U(X̂₀, Xₖ, ProjU, IndHₖ)

        ProjS_Xₖ = ProjS(Xₖ)
        iter += 1
        tired = iter > max_iter
        tolbestcirc = rel_error(Xₖ, ProjS_Xₖ)
        solved = tolbestcirc < tol 
        @info "iter = $iter, tolbestcirc = $tolbestcirc"
    end
    return Xₖ, iter
end



"""
    bestcirc_CRM(A::Matrix{T}, ProjS, ProjU) where {T<:Real}
    Using the Best Approximation circumcentered-reflection method to find the nearest correlation matrix to A
"""


function bestcirc_CRM2(X₀,
    ProjS,
    ProjU;
    max_iter::Int64=1000, tol::Float64=1e-7) where {T<:Real}
    
    X̂₀ = ProjU(X₀)
    ReflectU = x -> Reflection(x, ProjU)
    Xₖ = CRMiteration(X̂₀, ProjS(X̂₀), ReflectU)
    ProjS_Xₖ = ProjS(Xₖ)
    tired = false
    iter = 1
    tolbestcirc = norm(Xₖ - ProjS_Xₖ)
    solved = tolbestcirc < tol
    @info "iter = $iter, tolbestcirc = $tolbestcirc"
    while !(solved || tired)
        X̂ₖ = CRMiteration(Xₖ, ProjS_Xₖ, ReflectU)
        ProjS_X̂ₖ = ProjS(X̂ₖ)
        X̂ₖ_ProjS_X̂ₖ = X̂ₖ - ProjS_X̂ₖ
        X̂ₖ_ProjS_X̂ₖ /= norm(X̂ₖ_ProjS_X̂ₖ)
        IndHₖ = IndHalfspace(X̂ₖ_ProjS_X̂ₖ, dot(X̂ₖ_ProjS_X̂ₖ, ProjS_X̂ₖ))
        # Note that  Xₖ is the Projection of X̂₀ onto  Wₖ ∩ U
        Xₖ = Project_W_H_U(X̂₀, Xₖ, ProjU, IndHₖ)

        # end
        ProjS_Xₖ = ProjS(Xₖ)
        tolbestcirc = rel_error(Xₖ, ProjS_Xₖ)

        solved = tolbestcirc < tol
        iter += 1
        tired = iter > max_iter
        @info "iter = $iter, tolbestcirc = $tolbestcirc"

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
function Project_W_H_U(X₀, Xₖ, ProjU, IndH)
    T = eltype(X₀)
    IndW = IndHalfspace(X₀ - Xₖ, dot(X₀ - Xₖ, Xₖ))
    # Projection onto the intersection of U ∩ W is already Xₖ
    #  Check if Xₖ is in H. If so, it is the projection
    if IndH(Xₖ) == zero(T)
        return Xₖ
    end
    ProjectH = x -> ProjHyperplane(IndH.a, IndH.b, x)
    Project_U_cap_H = x -> Project_U_cap_Hyperplane(x, ProjU, ProjectH)
    Reflect_U_cap_H = x -> Reflection(x, Project_U_cap_H)
    ProjU_H_X₀ = Project_U_cap_H(X₀)
    # Return the circumcenter of X₀, R_W(X₀) and R_HR_W(X₀)
    Proj_U_W_H_X₀ = CRMiteration(X₀, Xₖ, Reflect_U_cap_H)
    vec_dist =  [X₀ - ProjU_H_X₀,
        X₀ - Proj_U_W_H_X₀]
    _, min_index = findmin(norm.(vec_dist))
    if min_index == 1 && IndW(ProjU_H_X₀) == zero(T)
        return ProjU_H_X₀
    else
        return Proj_U_W_H_X₀
    end
end

function Project_U_cap_Hyperplane(X, ProjU, ProjHyp)
    ReflectH = x -> Reflection(x, ProjHyp)
    ReflectU = x -> Reflection(x, ProjU)
    return CRMiteration(X, ReflectH, ReflectU)
end

function rel_error(a,b; p = 2)
    return norm(a - b, p)/norm(a, p)
end


##################################################################
### Examples
##################################################################


## Example 0 from Higham (2002)
A = Float64[1 1 0
            1 1 1
            0 1 1]
## Solving with Boyle-Dykstra's algorithm
@info "Solving with Boyle-Dykstra's algorithm"
X_BD, iter = dykstra_NCM(A, ProjS, ProjU, max_iter=1000, tol=1e-6)

norm(X_BD - A)

# q = [-0.4814, 0.7324, -0.4814]
# X_BD * q

## Solving with CRM
# @info "Solving with CRM"
# ResutlsCRM = CRM(A, ProjS, ProjU)
# X_CRM = ResutlsCRM.xApprox

# norm(X_CRM - A)


## Solving with bestcirc
@info "Solving with bestcirc_MAP"
X_BC_MAP, iter = bestcirc_MAP(A, ProjS, ProjU, tol=1e-6)

norm(X_BC_MAP  - A)


##Solving with bestcirc_Old
@info "Solving with bestcirc_CRM2"
X_BC_CRM, iter = bestcirc_CRM2(A, ProjS, ProjU, tol = 1e-6)

norm(X_BC_CRM - A)





## Example 1 from Higham (2002)
n = 2000
A = Tridiagonal(diagm(0 => fill(2.0,n), -1 => -ones(n-1), 1 => -ones(n-1)))
## Solving with Boyle-Dykstra's algorithm
@info "Solving with Boyle-Dykstra's algorithm"
X_BD, iter = dykstra_NCM(ProjU(A), ProjS, ProjU, tol = 1e-3)

norm(A - X_BD)

ResutlsCRM = CRM(A, ProjS, ProjU, EPSVAL = 1e-3)
X_CRM = ResutlsCRM.xApprox
norm(X_CRM - A) 



## Solving with bestcirc_MAP
@info "Solving with bestcirc_MAP"
X_BC_MAP, iter = bestcirc_MAP(A, ProjS, ProjU, tol = 1e-8)

norm(X_BC_MAP  - A)

## Solving with bestcirc_Old
@info "Solving with bestcirc_CRM2"
X_BC_CRM, iter = bestcirc_CRM2(A, ProjS, ProjU, tol = 1e-6)

norm(X_BC_CRM - A)