"""
    centralization!(xCRM, ProjectA, ProjectB)

Computes a centralized point
"""

function centralization!(xVEC::Vector,
    ProjectA::Function,
    ProjectB::Function)
    xVEC = MAP_iteration!(xVEC, ProjectA, ProjectB)
    xVEC = 0.5 * (xVEC + ProjectA(xVEC))
    return xVEC
end


"""
    centralization!(xCRM, ProjA, ProjectA, ProjectB)

Computes a centralized point
"""

function centralization!(xVEC::Vector,
    ProjA::Vector,
    ProjectA::Function,
    ProjectB::Function)
    xVEC = MAP_iteration!(xVEC, ProjA, ProjectB)
    xVEC = 0.5 * (xVEC + ProjectA(xVEC))
    return xVEC
end


"""
    parallelCRMiteration(xCRM, IsometryA, IsometryB)

Computes an iteration of the Cirumcentered-Reflection method
"""

function parallelCRMiteration!(xpCRM::Vector,
    IsometryA::Function,
    IsometryB::Function)
    xpCRM_RA = IsometryA(xpCRM)
    xpCRM_RB = IsometryB(xpCRM)
    xpCRM = parallelCRMiteration!(xpCRM, xpCRM_RA, xpCRM_RB)
    return xpCRM
end


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
    elseif xpCRM_RA ≈ xpCRM_RB
        xpCRM = FindCircumcentermSet([xpCRM, xpCRM_RA])
    else
        xpCRM = FindCircumcentermSet([xpCRM, xpCRM_RA, xpCRM_RB])
    end
    return xpCRM
end

"""
    centralizedCRM(x₀, ProjectA, ProjectB)

Centralized Cirumcentered-Reflection method
"""
function centralizedCRM(x₀::Vector, ProjectA::Function, ProjectB::Function;
    EPSVAL::Float64=1e-5, itmax::Int=100,
    filedir::String="", xSol::Vector=[],
    print_intermediate::Bool=false,
    gap_distance::Bool=true, isprod::Bool=false,
    centralized_onlyfirst::Bool=false)
    xcentCRM = x₀
    ReflectA(x) = Reflection(x, ProjectA)
    ReflectB(x) = Reflection(x, ProjectB)
    ProjA = ProjectA(xcentCRM)
    k = 0
    tolCRM = 1.0
    printOnFile(filedir, k, tolCRM, xcentCRM, deletefile=true, isprod=isprod)
    while tolCRM > EPSVAL && k < itmax
        print_intermediate ? printOnFile(filedir, 0, 0.0, ProjA, isprod=isprod) : nothing
        if gap_distance
            xcentCRM = centralization!(xcentCRM, ProjA, ProjectA, ProjectB)
            xcentCRM = parallelCRMiteration!(xcentCRM, ReflectA, ReflectB)
            ProjA = ProjectA(xcentCRM)
            ProjB = ProjectB(xcentCRM)
            tolCRM = max(norm(ProjA - xcentCRM), norm(ProjB - xcentCRM))
        else
            xcentCRMOld = copy(xcentCRM)
            xcentCRM = centralization!(xcentCRM, ProjA, ProjectA, ProjectB)
            xcentCRM = parallelCRMiteration!(xcentCRM, ReflectA, ReflectB)
            ProjA = ProjectA(xcentCRM)
            tolCRM = Tolerance(xcentCRM, xcentCRMOld, xSol)
        end
        k += 4
        printOnFile(filedir, k, tolCRM, xcentCRM, isprod=isprod)
    end
    isprod ? method = :CRMprod : method = :centCRM
    return Results(iter_total=k, final_tol=tolCRM, xApprox=xcentCRM, method=method)
end