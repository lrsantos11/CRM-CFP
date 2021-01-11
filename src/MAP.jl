"""
MAP
"""

function MAPiteration(xMAP::Vector, ProjectA, ProjectB)
    xMAP = ProjectA(xMAP)
    xMAP = ProjectB(xMAP)
    return xMAP  
end 

function MAP(x₀::Vector,ProjectA, ProjectB;EPSVAL::Float64=1e-5,itmax::Int = 100,filedir::String = "")
    k = 1
    tolMAP = 1.
    xMAP = x₀
    printoOnFile(filedir,xMAP',deletefile=true)
    while tolMAP > EPSVAL && k <= itmax
        xMAPOld = copy(xMAP)
        xMAP  = MAPiteration(xMAP, ProjectA, ProjectB)
        printoOnFile(filedir,xMAP')
        tolMAP = norm(xMAP-xMAPOld,Inf)
        k += 1
    end
    return xMAP, tolMAP, k
end