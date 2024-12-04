""" 
    Parallel Circumcentered Reflections Method
    
"""

"""
    parallelCRMiteration(xCRM, IsometryA, IsometryB)

Computes an iteration of the Cirumcentered-Reflection method
"""

function parallelCRMiteration!(xpCRM::Vector,
    IsometryA::Function,
    IsometryB::Function)
    xpCRM_RA = IsometryA(xpCRM)
    xpCRM_RB = IsometryB(xpCRM)
    xpCRM = parallelCRMiteration!(xpCRM, xpCRM_RA, xpCRM_RB)
    return xpCRM
end


"""
    parallelCRMiteration(xCRM, ReflectA, ReflectB)

Computes an iteration of the Cirumcentered-Reflection method
"""

function parallelCRMiteration!(xpCRM::Vector,
    xpCRM_RA::Vector,
    xpCRM_RB::Vector)
    if xpCRM_RA ≈ xpCRM
        xpCRM = FindCircumcentermSet([xpCRM, xpCRM_RB])
    elseif xpCRM_RB ≈ xpCRM
        xpCRM = FindCircumcentermSet([xpCRM, xpCRM_RA])
    elseif xpCRM_RA ≈ xpCRM_RB
        xpCRM = FindCircumcentermSet([xpCRM, xpCRM_RA])
    else
        xpCRM = FindCircumcentermSet([xpCRM, xpCRM_RA, xpCRM_RB])
    end
    return xpCRM
end


"""


    paralellCRMiteration!(xPCRM, Projections)

Computes a parallel CRM iteration
"""

using ThreadsX

function paralellCRMiteration!(xPCRM::Vector,
                        Projections)
    # Xaffine = [xPCRM]
    # for proj in Projections
        # push!(Xaffine, 2 * proj(xPCRM) - xPCRM)
    # end
    Xaffine = ThreadsX.map((proj) -> 2 * proj(xPCRM) - xPCRM, Projections)
    push!(Xaffine, xPCRM)
    return find_circumcenter!(xPCRM, Xaffine)
end




"""
    pCRM(x₀,Projections)

    Parallel Circumcentered Reflection Method
"""

function pCRM(x₀::AbstractVector,
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
    tolpCRM = 1.0
    xpCRM = copy(x₀)
    printOnFile(filedir, k, tolpCRM, xpCRM, deletefile=true, isprod=isprod)
    while tolpCRM > EPSVAL && k < itmax
        print_intermediate ? printOnFile(filedir, 0, 0.0, ProjA, isprod=isprod) : nothing
        if gap_distance
            @error "Not implemented yet"
        else
            xpCRMOld = copy(xpCRM)
            xpCRM = paralellCRMiteration!(xpCRM, Projections)
            tolpCRM = Tolerance(xpCRM, xpCRMOld, xSol)
        end
        k += num_sets
        printOnFile(filedir, k, tolpCRM, xpCRM, isprod=isprod)
    end
    isprod ? method = :pCRMprod : method = :pCRM
    return Results(iter_total=k, final_tol=tolpCRM, xApprox=xpCRM, method=method)

end