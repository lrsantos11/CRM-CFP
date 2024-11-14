""" 
    Circumcentered-Successive Reflections Method
    
    """


"""
    CSRM_iteration!(xCSRM, Projections)

Computes a SPM iteration
"""
function CSRM_iteration!(xCSRM::Vector,
                        Projections)
    Xaffine = map(x -> 2 * x - xCSRM, Projections(xCSRM))
    push!(Xaffine, xCSRM)
    return find_circumcenter!(xCSRM, Xaffine)
end




"""
    CSRM(x₀,Projections)

    Circumcentered Simultaneous Reflection Method
"""

function CSRM(x₀::AbstractVector,
    Projections,
    num_sets::Int;
    EPSVAL::Float64=1e-5,
    itmax::Int=100,
    filedir::String="",
    xSol::Vector=[],
    print_intermediate::Bool=false,
    gap_distance::Bool=false,
    isprod::Bool=false,
    kwargs...)
    k = 0
    tolCSRM = 1.0
    xCSRM = copy(x₀)
    printOnFile(filedir, k, tolCSRM, xCSRM, deletefile=true, isprod=isprod)
    while tolCSRM > EPSVAL && k < itmax
        print_intermediate ? printOnFile(filedir, 0, 0.0, ProjA, isprod=isprod) : nothing
        if gap_distance
            @error "Not implemented yet"
        else
            xCSRMOld = copy(xCSRM)
            xCSRM = CSRM_iteration!(xCSRM, Projections)
            tolCSRM = Tolerance(xCSRM, xCSRMOld, xSol)
        end
        k += num_sets
        printOnFile(filedir, k, tolCSRM, xCSRM, isprod=isprod)
    end
    isprod ? method = :SPMprod : method = :SPM
    return Results(iter_total=k, final_tol=tolCSRM, xApprox=xCSRM, method=method)

end