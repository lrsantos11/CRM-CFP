"""
    SPMiteration!(xSPM,ProjectA,ProjectB)

Computes a SPM iteration
"""
function SPM_iteration!(xSPM::Vector, 
                        ProjectionA::Function, 
                        ProjectionB::Function)
    xSPM = 0.5*(ProjectionA(xSPM) + ProjectionB(xSPM))                       
    return xSPM
end 

"""
    SPMiteration!(xSPM,ProjA,ProjectB)

Computes a SPM iteration
"""
function SPM_iteration!(xSPM::Vector, 
                        ProjA::Vector, 
                        ProjectionB::Function)
    xSPM = 0.5*(ProjA + ProjectionB(xSPM))                       
    return xSPM
end 



"""
    SPM(x₀,ProjectA, ProjectB)

    Simultaneous Projections Method
"""
function SPM(x₀::Vector, ProjectA::Function, ProjectB::Function; 
            EPSVAL::Float64 = 1e-5,
            itmax::Int = 100,
            filedir::String = "",
            xSol::Vector = [],
            print_intermediate::Bool = false,
            gap_distance::Bool = true, 
            isprod::Bool = false)
    k = 0
    tolSPM = 1.
    xSPM = x₀
    printOnFile(filedir,k, tolSPM, xSPM , deletefile=true, isprod=isprod)
    ProjA = ProjectA(xSPM)
    while tolSPM > EPSVAL && k < itmax
        print_intermediate ? printOnFile(filedir, 0, 0., ProjA , isprod=isprod) : nothing
        if gap_distance
            xSPM  = SPM_iteration!(xSPM, ProjA, ProjectB)
            ProjA = ProjectA(xSPM)
            tolSPM = norm(ProjA-xSPM)
        else
            xSPMOld = copy(xSPM)
            xSPM  = SPM_iteration!(xSPM, ProjA, ProjectB)
            ProjA = ProjectA(xSPM)
            tolSPM = Tolerance(xSPM,xSPMOld,xSol)
        end
        k += 2
        printOnFile(filedir,k, tolSPM, xSPM , isprod=isprod)
    end
    isprod ? method = :SPMprod : method = :SPM
    return Results(iter_total= k,final_tol=tolSPM,xApprox=xSPM,method=method)
end


