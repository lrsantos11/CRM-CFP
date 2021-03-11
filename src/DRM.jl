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
    EPSVAL::Float64=1e-5,itmax::Int = 100,filedir::String = "", xSol::Vector = [],
    print_intermediate::Bool=false,gap_distance::Bool=false)
    k = 0
    tolDRM = 1.
    xDRM = x₀
    ReflectA(x) = Reflection(x,ProjectA)
    ReflectB(x) = Reflection(x,ProjectB)
    printoOnFile(filedir,hcat(k, tolDRM, xDRM'),deletefile=true)
    while tolDRM > EPSVAL && k < itmax
        xDRMOld = copy(xDRM)
        print_intermediate ?  printoOnFile(filedir,hcat(nothing,nothing,ProjectA(xDRM)')) : nothing
        xDRM  = DRMiteration(xDRM, ReflectA, ReflectB)
        tolDRM = gap_distance ? norm(ProjectA(xDRM)-ProjectB(xDRM)) : Tolerance(xDRM,xDRMOld,xSol)
        k += 1
        printoOnFile(filedir,hcat(k, tolDRM, xDRM'))
    end
    return Results(iter_total= k,final_tol=tolDRM,xApprox=xDRM,method=:DRM)
end

"""
    DRMprod(x₀, SetsProjections)

Douglas–Rachford method on Pierra's product space reformulation
"""
function DRMprod(x₀::Vector{Float64},SetsProjections::Vector{Function}; 
    EPSVAL::Float64=1e-5,itmax::Int = 100,filedir::String = "", xSol::Vector = [],
    print_intermediate::Bool=false,gap_distance::Bool=false)
    k = 0
    tolDRMprod = 1.
    num_sets = length(SetsProjections)
    xDRMprod = Vector[]
    for i = 1:num_sets
        push!(xDRMprod,x₀)
    end
    ProjectAprod(x) = ProjectProdSets(x,SetsProjections)
    ProjectBprod(x) = ProjectProdDiagonal(x)
    ReflectA(x) = Reflection(x,ProjectAprod)
    ReflectB(x) = Reflection(x,ProjectBprod)
    printoOnFile(filedir,hcat(k, tolDRMprod, xDRMprod[1]'),deletefile=true)
    while tolDRMprod > EPSVAL && k < itmax
        xDRMprodOld = copy(xDRMprod)
        print_intermediate ?  printoOnFile(filedir,hcat(nothing,nothing,(ProjectA(xDRMprod))[1]')) : nothing
        xDRMprod  = DRMiteration(xDRMprod, ReflectA, ReflectB)
        tolDRMprod = gap_distance ? norm(ProjectAprod(xDRMprod)-xDRMprod) : Tolerance(xDRMprod,xDRMprodOld,xSol)
        k += 1
        printoOnFile(filedir,hcat(k, tolDRMprod, xDRMprod[1]'))
    end
    return Results(iter_total= k,
                  final_tol=tolDRMprod,xApprox=xDRMprod[1],method=:DRMprod)
end    


DRMprod(x₀::Vector,ProjectA::Function, ProjectB::Function; kwargs...) =   DRMprod(x₀,[ProjectA,ProjectB];kwargs...)

"""
    DRMprod(x₀, SetsIndicators)

Cirumcentered-Reflection method on Pierra's product space reformulation using
indicator functinos from `ProximalOperators.jl`
"""
function DRMprod(x₀::Vector{Float64},SetsIndicators::Vector{ProximableFunction}; 
    EPSVAL::Float64=1e-5,itmax::Int = 100,filedir::String = "", xSol::Vector = [],
    print_intermediate::Bool=false,gap_distance::Bool=false)
    k = 0
    tolDRMprod = 1.
    num_sets = length(SetsIndicators)
    xDRMprod = Vector[]
    for i = 1:num_sets
        push!(xDRMprod,x₀)
    end
    ProjectAprod(x) = ProjectSetsIndicators(x,SetsIndicators)
    ProjectBprod(x) = ProjectProdDiagonal(x)
    ReflectA(x) = Reflection(x,ProjectAprod)
    ReflectB(x) = Reflection(x,ProjectBprod)
    printoOnFile(filedir,hcat(k, tolDRMprod, xDRMprod[1]'),deletefile=true)
    while tolDRMprod > EPSVAL && k < itmax
        xDRMprodOld = copy(xDRMprod)
        print_intermediate ?  printoOnFile(filedir,hcat(nothing,nothing,(ProjectA(xDRMprod))[1]')) : nothing
        xDRMprod  = DRMiteration(xDRMprod, ReflectA, ReflectB)
        tolDRMprod = gap_distance ? norm(ProjectAprod(xDRMprod)-xDRMprod) : Tolerance(xDRMprod,xDRMprodOld,xSol)
        k += 1
        printoOnFile(filedir,hcat(k, tolDRMprod, xDRMprod[1]'))
    end
    return Results(iter_total= k,
                  final_tol=tolDRMprod,xApprox=xDRMprod[1],method=:DRMprod)
end 