"""
    CRMiteration(xCRM, ReflectA, ReflectB)

Computes an iteration of the Cirumcentered-Reflection method
"""

function CRMiteration(xCRM::Vector, 
                      ReflectA::Function, 
                      ReflectB::Function)
    xCRM_RA = ReflectA(xCRM)
    xCRM_RBRA = ReflectB(xCRM_RA)
    if xCRM_RA ≈ xCRM
        xCRM = FindCircumcentermSet([xCRM, xCRM_RBRA])
    elseif xCRM_RBRA ≈ xCRM_RA
        xCRM =FindCircumcentermSet([xCRM,  xCRM_RA])
    else
        xCRM = FindCircumcentermSet([xCRM, xCRM_RA, xCRM_RBRA])
    end
    return xCRM  
end 

"""
    CRMiteration(xCRM, ReflectA, ReflectB)

Computes an iteration of the Cirumcentered-Reflection method
"""

function CRMiteration(xCRM::Vector, 
                      ProjA::Vector, 
                      ReflectB::Function)
    xCRM_RA = 2*ProjA - xCRM
    xCRM_RBRA = ReflectB(xCRM_RA)
    if (xCRM_RA ≈ xCRM)
        xCRM = FindCircumcentermSet([xCRM, xCRM_RBRA])
    elseif (xCRM_RBRA ≈ xCRM_RA)
        xCRM = FindCircumcentermSet([xCRM,  xCRM_RA])
    else
        xCRM = FindCircumcentermSet([xCRM, xCRM_RA, xCRM_RBRA])
    end
    return xCRM  
end 

"""
    CRM(x₀, ProjectA, ProjectB)

Cirumcentered-Reflection method
"""
function CRM(x₀::Vector,ProjectA::Function, ProjectB::Function; 
             EPSVAL::Float64 = 1e-5,
             itmax::Int = 100,
             filedir::String = "", 
             xSol::Vector = [],
             print_intermediate::Bool = false,
             gap_distance::Bool = false,
             isprod::Bool = false, 
             verbose::Bool = false)
    xCRM = ProjectB(x₀)
    ReflectA(x) = Reflection(x,ProjectA)
    ReflectB(x) = Reflection(x,ProjectB)
    ProjA = ProjectA(xCRM)
    k = 0
    tolCRM = 1.
    printOnFile(filedir,k, tolCRM, xCRM ,deletefile=true, isprod=isprod)
    while tolCRM > EPSVAL && k < itmax
        print_intermediate ?  printOnFile(filedir,0, 0., ProjA,isprod=isprod) : nothing
        if gap_distance
            xCRM  = CRMiteration(xCRM, ProjA, ReflectB)
            ProjA = ProjectA(xCRM)                     
            ProjB = ProjectB(xCRM)
            tolCRM = max(norm(ProjA-xCRM),norm(ProjB-xCRM))       
        else
            xCRMOld = copy(xCRM)
            xCRM  = CRMiteration(xCRM, ProjA, ReflectB)
            ProjA = ProjectA(xCRM)
            tolCRM = Tolerance(xCRM,xCRMOld,xSol)
        end
        k += 2
        printOnFile(filedir,k, tolCRM, xCRM, isprod=isprod)
end
    isprod ? method = :CRMprod : method = :CRM
    return Results(iter_total = k,final_tol=tolCRM,xApprox=xCRM,method=method)
end

"""
    CRMprod(x₀, Projections)

Cirumcentered-Reflection method on Pierra's product space reformulation
"""
function CRMprod(x₀::Vector{Float64},Projections::Vector{Function}; kwargs...)
    xCRMprod = [x₀ for _ in Projections]
    ProjectAprod(x) = ProjectProdSpace(x,Projections)
    ProjectBprod(x) = ProjectProdDiagonal(x)
    results = CRM(xCRMprod, ProjectAprod, ProjectBprod, isprod=true; kwargs...)
    return results
end    

"""
    CRMprod(x₀, ProjectA, ProjectB)
"""
CRMprod(x₀::Vector{Float64},ProjectA::Function, ProjectB::Function;kwargs...) = CRMprod(x₀,[ProjectA,ProjectB];kwargs...) 