__precompile__()
include("auxfunctions.jl")
global const ZERO_VAL = 1e-15
include("CRM.jl")
include("MAP.jl")



####################################
"""
        ProjectEpiQuadratic(s,v,α)
Project ``(s,v) ∈ R^{n+1}`` onto the epigraph of ``f(x) = αx^Tx``, such that ``f(v) '≤ s``.
"""
function ProjectEpiQuadratic(s::Number,v::Union{AbstractArray,Number}; α::Float64=1.0)
        if α*dot(v,v) <= s
                return s, v
        end
        #PolynomialCoefficients
        # a3μ³ + a2μ² + a1μ  + a0 = 0    
        a0 = s-α*dot(v,v)
        a1 = 4*α*s + 1.
        a2 = 4*α^2*s + 4*α
        a3 = 4*α^2
        r = roots([a0,a1,a2,a3])
        indexreal =  findall(x->abs.(x)<1e-12,imag.(r))
        μ  = (maximum(real.(r[indexreal])))
        x = 1/(1+2*α*μ ) * v
        t =  μ + s
        return t, x
end
####################################
"""
        ProjectEpiQuadratic(x,α)
Project ``x = [x₀,t₀] ∈ R^{n+1}`` onto the epigraph of ``f(u) = αu^Tu``, such that ``f(x₀) '≤ t₀   ``.
"""
function ProjectEpiQuadratic(x::AbstractArray; α::Float64=1.0)
    t, x = ProjectEpiQuadratic(x[end],x[1:end-1], α=α)
    return [x;t]
end