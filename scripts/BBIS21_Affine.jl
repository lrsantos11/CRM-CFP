##################################################################
## Basic Functions for Linear Feasibility Problems  tests
##################################################################
##
"""
centralizedCRM(LPProblem)
Uses cCRM to find a point into intersection of two Affines 
"""
function centralizedCRM(mpsfile::String;
    kwargs...)
    nrow, ncol, cc, A, b, xlb, xub = LPDataFromMPS(mpsfile)
    Affine = IndAffine(A, b)
    Box = IndBox(xlb, xub)
    Random.seed!(1234)
    x₀ = StartingPoint(ncol)
    ProjectA(x) = ProjectIndicator(Box, x)
    ProjectB(x) = ProjectIndicator(Affine, x)
    @info "cCRM for LP"
    return centralizedCRM(x₀, ProjectA, ProjectB; kwargs...)
end

"""

MAP(x₀, Affines)
Uses MAP to find  a point into intersection of two Affines 
"""
function MAP(mpsfile::String;
    kwargs...)
    nrow, ncol, cc, A, b, xlb, xub = LPDataFromMPS(mpsfile)
    Affine = IndAffine(A, b)
    Box = IndBox(xlb, xub)
    Random.seed!(1234)
    x₀ = StartingPoint(ncol)
    ProjectA(x) = ProjectIndicator(Box, x)
    ProjectB(x) = ProjectIndicator(Affine, x)
    @info "MAP for LP"
    return MAP(x₀, ProjectA, ProjectB; kwargs...)
end


"""
SPM(x₀, Affines)
Uses SPM to find  a point into intersection of two Affines 
"""
function SPM(mpsfile::String;
    kwargs...)
    nrow, ncol, cc, A, b, xlb, xub = LPDataFromMPS(mpsfile)
    Affine = IndAffine(A, b)
    Box = IndBox(xlb, xub)
    Random.seed!(1234)
    x₀ = StartingPoint(ncol)
    ProjectA(x) = ProjectIndicator(Box, x)
    ProjectB(x) = ProjectIndicator(Affine, x)
    @info "SPM for LP"
    return SPM(x₀, ProjectA, ProjectB; kwargs...)
end
##

"""
CRM(x₀, Affines)
Uses CRM to find  a point into intersection of two Affines 
"""
function CRM(mpsfile::String;
    kwargs...)
    nrow, ncol, cc, A, b, xlb, xub = LPDataFromMPS(mpsfile)
    Affine = IndAffine(A, b)
    Box = IndBox(xlb, xub)
    Random.seed!(1234)
    x₀ = StartingPoint(ncol)
    ProjectA(x) = ProjectIndicator(Box, x)
    ProjectB(x) = ProjectIndicator(Affine, x)
    @info "CRM for LP"
    return CRM(x₀, ProjectA, ProjectB; kwargs...)
end

##




# mpsfile = NETLIB_PRE_DIR * "/afiro"
mpsfile = NETLIB_DIR * "/AFIRO.sif"
itmax = 400
ε = 1e-6
##
xCRM = CRM(mpsfile, gap_distance=true, EPSVAL=ε, itmax=itmax)
##
xcCRM = centralizedCRM(mpsfile, gap_distance=true, EPSVAL=ε, itmax=itmax)
##
xMAP = MAP(mpsfile, gap_distance=true, EPSVAL=ε, itmax=itmax)
