"""
    MAPiteration(xMAP,ProjectA,ProjectB)

Computes a MAP iteration
"""

function MAPiteration(xMAP::Vector, ProjectA::Function, ProjectB::Function,
                    filedir::String="",print_intermediate::Bool=true)
    xMAP = ProjectA(xMAP)
    print_intermediate ? printoOnFile(filedir,hcat(nothing, nothing, xMAP')) : nothing
    xMAP = ProjectB(xMAP)
    return xMAP  
end 


"""
    MAP(x₀,ProjectA, ProjectB)

    Method of Alternating Projections
"""
function MAP(x₀::Vector,ProjectA::Function, ProjectB::Function; 
        EPSVAL::Float64=1e-5,itmax::Int = 100,filedir::String = "",xSol::Vector = [],
        print_intermediate::Bool=true,gap_distance::Bool=false)
    k = 0
    tolMAP = 1.
    xMAP = x₀
    printoOnFile(filedir,hcat(k, tolMAP, xMAP'),deletefile=true)
    while tolMAP > EPSVAL && k < itmax
        xMAPOld = copy(xMAP)
        xMAP  = MAPiteration(xMAP, ProjectA, ProjectB,filedir,print_intermediate)
        tolMAP = gap_distance ? norm(ProjectA(xMAP)-xMAP) : Tolerance(xMAP,xMAPOld,xSol)
        k += 1
        printoOnFile(filedir,hcat(k, tolMAP, xMAP'))
    end
    return Results(iter_total= k,final_tol=tolMAP,xApprox=xMAP,method=:MAP)
end


"""
    MAPprod(x₀, SetsProjections)

    Method of Alternating projections on Pierra's product space reformulation
"""
function MAPprod(x₀::Vector{Float64},Projections::Vector; 
    EPSVAL::Float64=1e-5,itmax::Int = 100,filedir::String = "", xSol::Vector = [],
    print_intermediate::Bool=false,gap_distance::Bool=false)
    k = 0
    tolMAPprod = 1.
    num_sets = length(Projections)
    xMAPprod = Vector[]
    for i = 1:num_sets
        push!(xMAPprod,x₀)
    end
    ProjectAprod(x) = ProjectProdSpace(x,Projections)
    ProjectBprod(x) = ProjectProdDiagonal(x)
    printoOnFile(filedir,hcat(k, tolMAPprod, xMAPprod[1]'),deletefile=true)
    while tolMAPprod > EPSVAL && k < itmax
        xMAPprodOld = copy(xMAPprod)
        xMAPprod  = MAPiteration(xMAPprod, ProjectAprod, ProjectBprod,filedir,print_intermediate)
        tolMAPprod = gap_distance ? norm(ProjectAprod(xMAPprod)-xMAPprod) : Tolerance(xMAPprod,xMAPprodOld,xSol)
        k += 1
        printoOnFile(filedir,hcat(k, tolMAPprod, xMAPprod[1]'))
    end
    return Results(iter_total= k,
                  final_tol=tolMAPprod,xApprox=xMAPprod[1],method=:MAPprod)
end    

"""
    MAPprod(x₀, ProjectA, ProjectB)
"""
MAPprod(x₀::Vector{Float64},ProjectA::Function, ProjectB::Function;kwargs...) = MAPprod(x₀,[ProjectA,ProjectB],kwargs...) 