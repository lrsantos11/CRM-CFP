"""
    CRMiteration(xCRM, ReflectA, ReflectB)

Computes an iteration of the Cirumcentered-Reflection method
"""

function CRMiteration(xCRM::Vector, ReflectA::Function, ReflectB::Function)
    xCRM_RA = ReflectA(xCRM)
    xCRM_RBRA = ReflectB(xCRM_RA)
    if norm(xCRM_RA - xCRM)<ZERO_VAL
        xCRM = FindCircumcentermSet([xCRM, xCRM_RBRA])
    elseif norm(xCRM_RBRA - xCRM_RA)<ZERO_VAL
        xCRM =FindCircumcentermSet([xCRM,  xCRM_RA])
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
    EPSVAL::Float64=1e-5,itmax::Int = 100,filedir::String = "", xSol::Vector = [],print_intermediate::Bool=false)
    k = 0
    tolCRM = 1.
    xCRM = x₀
    ReflectA(x) = Reflection(x,ProjectA)
    ReflectB(x) = Reflection(x,ProjectB)
    printoOnFile(filedir,xCRM',deletefile=true)
    while tolCRM > EPSVAL && k < itmax
        xCRMOld = copy(xCRM)
        print_intermediate ?  printoOnFile(filedir,ProjectA(xCRM)') : nothing
        xCRM  = CRMiteration(xCRM, ReflectA, ReflectB)
        printoOnFile(filedir,xCRM')
        tolCRM = Tolerance(xCRM,xCRMOld,xSol)
        k += 1
    end
    return Results(iter_total= k,final_tol=tolCRM,xApprox=xCRM,method="CRM")
end
"""
    CRMprod(x₀, SetsProjections)

Cirumcentered-Reflection method on Pierra's product space reformulation
"""
function CRMprod(x₀::Vector{Float64},SetsProjections::Vector{Function}; 
    EPSVAL::Float64=1e-5,itmax::Int = 100,filedir::String = "", xSol::Vector = [],print_intermediate::Bool=false)
    k = 0
    tolCRMprod = 1.
    # dim_space = length(x₀)
    num_sets = length(SetsProjections)
    xCRMprod = Vector[]
    # sizehint!(xCRMprod,num_sets*dim_space)
    for i = 1:num_sets
        push!(xCRMprod,x₀)
    end
    ProjectAprod(x) = ProjectProdSets(x,SetsProjections)
    ProjectBprod(x) = ProjectProdDiagonal(x,num_sets)
    ReflectA(x) = Reflection(x,ProjectAprod)
    ReflectB(x) = Reflection(x,ProjectBprod)
    printoOnFile(filedir,xCRMprod[1]',deletefile=true)
    while tolCRMprod > EPSVAL && k < itmax
        xCRMprodOld = copy(xCRMprod)
        print_intermediate ?  printoOnFile(filedir,(ProjectA(xCRMprod))[1]') : nothing
        xCRMprod  = CRMiteration(xCRMprod, ReflectA, ReflectB)
        printoOnFile(filedir,xCRMprod[1]')
        tolCRMprod = Tolerance(xCRMprod,xCRMprodOld,xSol)
        k += 1
    end
    return Results(iter_total= k,
                  final_tol=tolCRMprod,xApprox=xCRMprod[1],method="CRMprod")
end    


CRMprod(x₀::Vector,ProjectA::Function, ProjectB::Function; kwargs...) =   CRMprod(x₀,[ProjectA,ProjectB];kwargs...)