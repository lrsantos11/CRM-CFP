

"""
    parallelCRMiteration(xCRM, ReflectA, ReflectB)

Computes an iteration of the Cirumcentered-Reflection method
"""

function parallelCRMiteration!(xpCRM::Vector,
    xpCRM_RA::Vector,
    xpCRM_RB::Vector)
    if xpCRM_RA ≈ xpCRM
        xpCRM = FindCircumcentermSet([xpCRM, xpCRM_RB])
    elseif xpCRM_RB ≈ xpCRM
        xpCRM = FindCircumcentermSet([xpCRM, xpCRM_RA])
    else
        xpCRM = FindCircumcentermSet([xpCRM, xpCRM_RA, xpCRM_RB])
    end
    return xpCRM
end




"""
    SucCentCRM_Cyclic(x₀, Projections)

Successive Centralized Cirumcentered-Reflection method using Cyclic control sequence
"""
function SucCentCRM_Cyclic(x₀::Vector,
    Projections::Vector{Function};
    EPSVAL::Float64=1e-6, itmax::Int=100,
    filedir::String="",
    verbose::Bool=false,
    error_sum::Bool=true
)
    xSucCentCRM = copy(x₀)
    m = length(Projections)
    k = 0
    tol = 1.0
    printOnFile(filedir, k, tol, xSucCentCRM, deletefile=true)
    solved = false
    tired = false
    while !(solved || tired)
        xSucCentCRMOld = copy(xSucCentCRM)
        @inbounds for i in 1:m
            ProjFirst = Projections[i](xSucCentCRM)
            iplus1 = i == m ? 1 : i + 1 #  moves to the first one if it is the last set
            xSucCentCRM = centralization!(xSucCentCRM, ProjFirst, Projections[i], Projections[iplus1])
            ReflectA = 2 * Projections[i](xSucCentCRM) - xSucCentCRM
            ReflectB = 2 * Projections[iplus1](xSucCentCRM) - xSucCentCRM
            xSucCentCRM = parallelCRMiteration!(xSucCentCRM, ReflectA, ReflectB)
        end
        k += 1
        solved, tol = stopping_criteria(xSucCentCRM, xSucCentCRMOld, Projections, EPSVAL, error_sum, verbose)
        tired = k ≥ itmax
        printOnFile(filedir, k, tol, xSucCentCRM)
    end
    return Results(iter_total=k, proj_total = 4*m*k,  final_tol=tol, xApprox=xSucCentCRM, method=:SucCentCRM_Cyclic)
end




"""
    SucCentCRM_AlmostViolatedDist(x₀, Projections)

Successive Centralized Cirumcentered-Reflection method using Cyclic control sequence distance version
"""
function SucCentCRM_AlmostViolatedDist(x₀::Vector,
    Projections::Vector{Function};
    EPSVAL::Float64=1e-6, itmax::Int=100,
    filedir::String="",
    verbose::Bool=false,
    error_sum::Bool=true
)
    xSucCentCRMViolDist = copy(x₀)
    k = 0
    tol = 1.0
    printOnFile(filedir, k, tol, xSucCentCRMViolDist, deletefile=true)
    solved = false
    tired = false
    while !(solved || tired)
        xSucCentCRMViolDistOld = copy(xSucCentCRMViolDist)
        Projℓ = [proj(xSucCentCRMViolDist) for proj ∈ Projections]
        _, ℓₖ = findmax(norm.(Ref(x₀) .- Projℓ))
        Projℓ_z = @views Projℓ[ℓₖ]
        Projr = [proj(Projℓ_z) for proj ∈ Projections]
        _, rₖ = findmax(norm.(Ref(Projℓ_z) .- Projr))
        xSucCentCRMViolDist = centralization!(xSucCentCRMViolDist, Projℓ_z, Projections[ℓₖ], Projections[rₖ])
        ReflectA = 2 * Projections[ℓₖ](xSucCentCRMViolDist) - xSucCentCRMViolDist
        ReflectB = 2 * Projections[rₖ](xSucCentCRMViolDist) - xSucCentCRMViolDist
        xSucCentCRMViolDist = parallelCRMiteration!(xSucCentCRMViolDist, ReflectA, ReflectB)
        k += 1
       solved, tol = stopping_criteria(xSucCentCRMViolDist,  xSucCentCRMViolDistOld, Projections, EPSVAL, error_sum, verbose)
       tired = k ≥ itmax
        printOnFile(filedir, k, tol, xSucCentCRMViolDist)
    end
    return Results(iter_total=k, proj_total = k*(2m+4),final_tol=tol, xApprox=xSucCentCRMViolDist, method=:xSucCentCRMViolDist)
end

"""
    SucCentCRM_AlmostViolatedFunc(x₀, Projections)

Successive Centralized Cirumcentered-Reflection method using Cyclic control sequence functional version
"""
function SucCentCRM_AlmostViolatedFunc(x₀::Vector,
    Projections::Vector{Function},
    FunctionEval::Vector{Function};   
    EPSVAL::Float64=1e-6, itmax::Int=100,
    filedir::String="",
    verbose::Bool=false,
    error_sum::Bool=true
)
    xSucCentCRMViolFunc = copy(x₀)
    k = 0
    tol = 1.0
    printOnFile(filedir, k, tol, xSucCentCRMViolFunc, deletefile=true)
    solved = false
    tired = false
    while !(solved || tired)
        xSucCentCRMViolFuncOld = copy(xSucCentCRMViolFunc)
        _, ℓₖ = findmax([func(xSucCentCRMViolFunc) for func ∈ FunctionEval])
        Projℓ_z = Projections[ℓₖ](xSucCentCRMViolFunc)
        _, rₖ = findmax([func(Projℓ_z) for func ∈ FunctionEval])
        xSucCentCRMViolFunc = centralization!(xSucCentCRMViolFunc, Projℓ_z, Projections[ℓₖ], Projections[rₖ])
        ReflectA = 2 * Projections[ℓₖ](xSucCentCRMViolFunc) - xSucCentCRMViolFunc
        ReflectB = 2 * Projections[rₖ](xSucCentCRMViolFunc) - xSucCentCRMViolFunc
        xSucCentCRMViolFunc = parallelCRMiteration!(xSucCentCRMViolFunc, ReflectA, ReflectB)
        k += 1
        solved, tol = stopping_criteria(xSucCentCRMViolFunc, xSucCentCRMViolFuncOld, Projections, EPSVAL, error_sum, verbose)
        tired = k ≥ itmax
        printOnFile(filedir, k, tol, xSucCentCRMViolFunc)
    end
    m = length(Projections)
    return Results(iter_total=k, proj_total = 5*k, final_tol=tol, xApprox=xSucCentCRMViolFunc, method=:xSucCentCRMViolFunc)
end



"""
    SePM(x₀, Projections)

Sequential Projections  Method 
"""
function SePM(x₀::Vector,
    Projections::Vector{Function};
    EPSVAL::Float64=1e-6, itmax::Int=100,
    filedir::String="",
    verbose::Bool=false,
    error_sum::Bool=true
)
    xSePM = copy(x₀)
    m = length(Projections)
    k = 0
    tol = 1.0
    printOnFile(filedir, k, tol, xSePM, deletefile=true)
    solved = false
    tired = false
    while !(solved || tired)
        xSePMOld = copy(xSePM)
        for i in 1:m
            @inbounds xSePM = Projections[i](xSePM)
        end
        k += 1
        solved, tol = stopping_criteria(xSePM, xSePMOld, Projections, EPSVAL, error_sum, verbose)
        tired = k ≥ itmax
        printOnFile(filedir, k, tol, xSePM)
    end
    return Results(iter_total=k, proj_total = k*m,final_tol=tol, xApprox=xSePM, method=:SePM)
end


"""
    stopping_criteria(x, x_old, EPSVAL, itmax)
"""
function stopping_criteria(x, x_old, Projections, EPSVAL, error_sum, verbose)
        current_error = error_sum ? sum([norm(proj(x) - x) for proj in Projections]) : 0.0
        tol = max(current_error, norm(x - x_old))
        verbose && @info "tol = $tol"
        solved = tol < EPSVAL
    return solved, tol
end
