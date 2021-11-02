"""
    centralization!(xCRM, ProjectA, ProjectB)

Computes a centralized point
"""

function centralization!(xVEC::Vector,
                        ProjectA::Function,
                        ProjectB::Function)
    xVEC = MAP_iteration!(xVEC,ProjectA,ProjectB)
    xVEC  = 0.5*(xVEC + ProjectA(xVEC))
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
    xVEC = MAP_iteration!(xVEC,ProjA,ProjectB)
    xVEC  = 0.5*(xVEC + ProjectA(xVEC))
    return xVEC
end


"""
    parallelCRMiteration(xCRM, IsometryA, IsometryB)

Computes an iteration of the Cirumcentered-Reflection method
"""

function parallelCRMiteration!(xCRM::Vector, 
                               IsometryA::Function, 
                               IsometryB::Function)
    xCRM_RA = IsometryA(xCRM)
    xCRM_RB = IsometryB(xCRM)
    if xCRM_RA ≈ xCRM
        xCRM = FindCircumcentermSet([xCRM, xCRM_RB])
    elseif xCRM_RB ≈ xCRM_RA
        xCRM = FindCircumcentermSet([xCRM,  xCRM_RA])
    else
        xCRM = FindCircumcentermSet([xCRM, xCRM_RA, xCRM_RB])
    end
    return xCRM  
end 

"""
    centralizedCRM(x₀, ProjectA, ProjectB)

Centralized Cirumcentered-Reflection method
"""
function centralizedCRM(x₀::Vector,ProjectA::Function, ProjectB::Function; 
                        EPSVAL::Float64=1e-5, itmax::Int = 100,
                        filedir::String = "", xSol::Vector = [],
                        print_intermediate::Bool = false,
                        gap_distance::Bool = true, isprod :: Bool = false,
                        centralized_onlyfirst :: Bool = false)
    xCRM = x₀
    ReflectA(x) = Reflection(x,ProjectA)
    ReflectB(x) = Reflection(x,ProjectB)
    ProjA = ProjectA(xCRM)
    k = 0
    tolCRM = 1.
    printOnFile(filedir, k, tolCRM, xCRM , deletefile=true, isprod=isprod)
    while tolCRM > EPSVAL && k < itmax
        print_intermediate ?  printOnFile(filedir,0, 0., ProjA,isprod=isprod) : nothing
        if gap_distance
            xCRM = centralization!(xCRM, ProjA, ProjectA, ProjectB)
            xCRM = parallelCRMiteration!(xCRM, ReflectA, ReflectB)
            ProjA = ProjectA(xCRM)
            tolCRM = norm(ProjA-xCRM)
        else
            xCRMOld = copy(xCRM)
            xCRM = centralization!(xCRM, ProjA, ProjectA, ProjectB)
            xCRM = parallelCRMiteration!(xCRM, ReflectA, ReflectB)
            ProjA = ProjectA(xCRM)
            tolCRM = Tolerance(xCRM,xCRMOld,xSol)
        end
        k += 4
        printOnFile(filedir,k, tolCRM, xCRM, isprod=isprod)
    end
    isprod ? method = :CRMprod : method = :centCRM
    return Results(iter_total= k,final_tol=tolCRM,xApprox=xCRM,method=method)
end