"""
DRM iteration
"""

function DRMiteration(xDRM::Vector, ReflectA::Function, ReflectB::Function)
    xDRM_RA = ReflectA(xDRM)
    xDRM_RBRA = ReflectB(xDRM_RA)
    xDRM = 0.5*(xDRM + xDRM_RBRA)
    return xDRM  
end 

"""
    DRM(x₀, ProjectA, ProjectB)

Douglas–Rachford method on Pierra's product space reformulation
"""

function DRM(x₀::Vector,ProjectA::Function, ProjectB::Function; 
    EPSVAL::Float64=1e-5,itmax::Int = 100,filedir::String = "", xSol::Vector = [],
    print_intermediate::Bool=false,gap_distance::Bool=false)
    k = 0
    tolDRM = 1.
    xDRM = x₀
    ReflectA(x) = Reflection(x,ProjectA)
    ReflectB(x) = Reflection(x,ProjectB)
    printOnFile(filedir,hcat(k, tolDRM, xDRM'),deletefile=true)
    while tolDRM > EPSVAL && k < itmax
        xDRMOld = copy(xDRM)
        print_intermediate ?  printOnFile(filedir,hcat(nothing,nothing,ProjectA(xDRM)')) : nothing
        xDRM  = DRMiteration(xDRM, ReflectA, ReflectB)
        tolDRM = gap_distance ? norm(ProjectA(xDRM)-ProjectB(xDRM)) : Tolerance(ProjectA(xDRM),ProjectA(xDRMOld),xSol)
        k += 2
        printOnFile(filedir,hcat(k, tolDRM, xDRM'))
    end
    return Results(iter_total= k,final_tol=tolDRM,xApprox=xDRM,method=:DRM)
end

"""
    DRMprod(x₀, Projections)

Douglas–Rachford method on Pierra's product space reformulation
"""
function DRMprod(x₀::Vector{Float64},Projections::Vector; 
    EPSVAL::Float64=1e-5,itmax::Int = 100,filedir::String = "", xSol::Vector = [],
    print_intermediate::Bool=false,gap_distance::Bool=false)
    k = 0
    tolDRMprod = 1.
    num_sets = length(Projections)
    xDRMprod = Vector[]
    for i = 1:num_sets
        push!(xDRMprod,x₀)
    end
    ProjectAprod(x) = ProjectProdSpace(x,Projections)
    ProjectBprod(x) = ProjectProdDiagonal(x)
    ReflectA(x) = Reflection(x,ProjectAprod)
    ReflectB(x) = Reflection(x,ProjectBprod)
    printOnFile(filedir,hcat(k, tolDRMprod, xDRMprod[1]'),deletefile=true)
    while tolDRMprod > EPSVAL && k < itmax
        xDRMprodOld = copy(xDRMprod)
        print_intermediate ?  printOnFile(filedir,hcat(nothing,nothing,(ProjectAprod(xDRMprod))[1]')) : nothing
        xDRMprod  = DRMiteration(xDRMprod, ReflectA, ReflectB)
        tolDRMprod = gap_distance ? norm(ProjectAprod(xDRMprod)-ProjectBprod(xDRMprod)) : Tolerance(xDRMprod,xDRMprodOld,xSol)
        k += 1
        printOnFile(filedir,hcat(k, tolDRMprod, xDRMprod[1]'))
    end
    return Results(iter_total= k,
                  final_tol=tolDRMprod,xApprox=xDRMprod[1],method=:DRMprod)
end    


"""
    DRMprod(x₀, ProjectA, ProjectB)
"""
DRMprod(x₀::Vector{Float64},ProjectA::Function, ProjectB::Function;kwargs...) = DRMprod(x₀,[ProjectA,ProjectB],kwargs...) 