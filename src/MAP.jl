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