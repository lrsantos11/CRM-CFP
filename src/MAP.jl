"""
MAP
"""

function MAPiteration(xMAP::Vector, ProjectA, ProjectB,filedir,print_intermediate)
    xMAP = ProjectA(xMAP)
    print_intermediate ? printoOnFile(filedir,xMAP') : nothing
    xMAP = ProjectB(xMAP)
    printoOnFile(filedir,xMAP')
    return xMAP  
end 

function MAP(x₀::Vector,ProjectA::Function, ProjectB::Function; 
        EPSVAL::Float64=1e-5,itmax::Int = 100,filedir::String = "",xSol::Vector = [],print_intermediate::Bool=true)
    k = 0
    tolMAP = 1.
    xMAP = x₀
    printoOnFile(filedir,xMAP',deletefile=true)
    while tolMAP > EPSVAL && k <= itmax
        xMAPOld = copy(xMAP)
        xMAP  = MAPiteration(xMAP, ProjectA, ProjectB,filedir,print_intermediate)
        tolMAP = Tolerance(xMAP,xMAPOld,xSol)
        k += 1
    end
    return Results(iter_total= k,final_tol=tolMAP,xApprox=xMAP,method="MAP")
end