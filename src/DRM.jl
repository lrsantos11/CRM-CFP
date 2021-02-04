"""
DRM
"""

function DRMiteration(xDRM::Vector, ReflectA::Function, ReflectB::Function)
    xDRM_RA = ReflectA(xDRM)
    xDRM_RBRA = ReflectB(xDRM_RA)
    xDRM = 0.5*(xDRM + xDRM_RBRA)
    return xDRM  
end 


function DRM(x₀::Vector,ProjectA::Function, ProjectB::Function; 
    EPSVAL::Float64=1e-5,itmax::Int = 100,filedir::String = "", xSol::Vector = [],print_intermediate::Bool=false)
    k = 0
    tolDRM = 1.
    xDRM = x₀
    ReflectA(x) = Reflection(x,ProjectA)
    ReflectB(x) = Reflection(x,ProjectB)
    printoOnFile(filedir,xDRM',deletefile=true)
    while tolDRM > EPSVAL && k < itmax
        xDRMOld = copy(xDRM)
        print_intermediate ?  printoOnFile(filedir,ProjectA(xDRM)') : nothing
        xDRM  = DRMiteration(xDRM, ReflectA, ReflectB)
        printoOnFile(filedir,xDRM')
        tolDRM = Tolerance(xDRM,xDRMOld,xSol)
        k += 1
    end
    return Results(iter_total= k,final_tol=tolDRM,xApprox=xDRM,method="DRM")
end