"""
MAP
"""

function MAPiteration(xMAP::Vector, ProjectA, ProjectB,filedir)
    xMAP = ProjectA(xMAP)
    printoOnFile(filedir,xMAP')
    xMAP = ProjectB(xMAP)
    printoOnFile(filedir,xMAP')
    return xMAP  
end 

function MAP(x₀::Vector,ProjectA::Function, ProjectB::Function; EPSVAL::Float64=1e-5,itmax::Int = 100,filedir::String = "",xSol::Vector = [])
    k = 1
    tolMAP = 1.
    xMAP = x₀
    printoOnFile(filedir,xMAP',deletefile=true)
    while tolMAP > EPSVAL && k <= itmax
        xMAPOld = copy(xMAP)
        xMAP  = MAPiteration(xMAP, ProjectA, ProjectB,filedir)
        tolMAP = Tolerance(xMAP,xMAPOld,xSol)
        k += 1
    end
    return Results(iter_total= k,final_tol=tolMAP,xApprox=xMAP,method="MAP")
end