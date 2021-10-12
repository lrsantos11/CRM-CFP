##################################################################
## Basic Functions for Two AffineSubspaces tests
##################################################################
const VectorAffine = Vector{ProximalOperators.IndAffineDirect{LinearAlgebra.QRCompactWY{Float64, Matrix{Float64}}, Float64, Float64, Matrix{Float64}, Vector{Float64}}} 
"""
AffineRandomPair(n,cmax,affine)

Generates a pair of random affine spaces with intersection 
of dimension `cmax`.
"""
function AffineRandomPair(n::Int,
                            cmax::Int;
                            affine::Bool = false)
        cmax >=2 ||  error("cmax must be greater than 3")
        n >= cmax ||  error("n must be greater than intersection dimension")
        @show mcex = rand(1:cmax) ## number of common normals
        maex = rand(1:n-3*mcex) ## number of extra normals of A
        mbex = rand(1:n-2*mcex-maex) ## number of extra normals of B

        Cex = randn(mcex,n)
        while !(rank(Cex) == mcex)
            Cex = randn(mcex,n)
        end

        Aex = randn(maex,n)
        while !(rank([Aex; Cex]) == maex +  mcex)
            Aex = randn(maex,n)
        end

        Bex = randn(mbex,n)
        while !(rank([Aex; Cex; Bex]) == maex +  mcex + mbex)
            Bex = randn(mbex,n)
        end

        ma = maex + mcex
        mb = mbex + mcex

        # Space A
        A = [Aex; Cex]

        # Space B
        B = [Cex; Bex]

        a = zeros(ma)
        b = zeros(mb)

        AffineA = IndAffine(A,a)
        AffineB = IndAffine(B,b)
        Affines = [AffineA, AffineB]

        return Affines
    end
"""
centralizedCRM(x₀, Affines)
Uses cCRM to find a point into intersection of two Affines 
"""
function centralizedCRM(x₀::Vector,
                        Affines::VectorAffine; 
                        kwargs...)
    ProjectA(x) = ProjectIndicator(Affines[1], x)
    ProjectB(x) = ProjectIndicator(Affines[2], x)
    return centralizedCRM(x₀, ProjectA, ProjectB; kwargs...)
end

"""

MAP(x₀, Affines)
Uses MAP to find  a point into intersection of two Affines 
"""
function MAP(x₀::Vector,
             Affines::VectorAffine; 
             kwargs...)
    ProjectA(x) = ProjectIndicator(Affines[1], x)
    ProjectB(x) = ProjectIndicator(Affines[2], x)
    return MAP(x₀, ProjectA, ProjectB; kwargs...)
end


"""
SPM(x₀, Affines)
Uses SPM to find  a point into intersection of two Affines 
"""
function SPM(x₀::Vector,
             Affines::VectorAffine; 
             kwargs...)
    ProjectA(x) = ProjectIndicator(Affines[1], x)
    ProjectB(x) = ProjectIndicator(Affines[2], x)
    return SPM(x₀, ProjectA, ProjectB; kwargs...)
end


"""
CRM(x₀, Affines)
Uses CRM to find  a point into intersection of two Affines 
"""
function CRM(x₀::Vector,
             Affines::VectorAffine; 
             kwargs...)
    ProjectA(x) = ProjectIndicator(Affines[1],x)
    ProjectB(x) = ProjectIndicator(Affines[2],x)
    x₀ =  ProjectB(x₀)
    return CRM(x₀, ProjectA, ProjectB; kwargs...)
end


"""
DRM(x₀, Affines)
Uses CRM to find  a point into intersection of two Affines 
"""
function DRM(x₀::Vector,
             Affines::VectorAffine; 
             kwargs...)
    ProjectA(x) = ProjectIndicator(Affines[1],x)
    ProjectB(x) = ProjectIndicator(Affines[2],x)
    x₀ =  ProjectB(x₀)
    return DRM(x₀, ProjectA, ProjectB; kwargs...)
end



##

n = 100
cmax = rand(1:ceil(Integer,n/5))

Affines = AffineRandomPair(n,cmax)

x₀ = StartingPoint(n)


@show CRM(x₀, Affines, gap_distance = true)
@show centralizedCRM(x₀, Affines, gap_distance = true);
