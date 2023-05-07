##################################################################
### Creating Ellipsoids Samples
##################################################################

"""
SampleTwoEllipsoids(n, p; λ=1.1, γ=1.5)
Creates two ellipsoids in ℝⁿ that intersect. The intersection is regulated by λ.

"""
function SampleTwoEllipsoids(n::Int,  # dimension
    p::Real; # sparsity of matrix A
    λ::Real=1.1, # parameter for touching ellipsoid
    γ::Real=1.5 # parameter for making A positive definite
)
    Ellipsoids = EllipsoidCRM[]
    A = Matrix(sprandn(n, n, p))
    A = (γ * I + A' * A)
    a = rand(n)
    b = A * a
    adotAa = dot(a, b)
    b .*= -1.0
    α = (1 + γ) * adotAa
    push!(Ellipsoids, EllipsoidCRM(A, b, α))
    TouchEll, Center2, TouchPoint = GenerateTouchingEllipsoid(Ellipsoids[1], n, λ=λ)
    push!(Ellipsoids, TouchEll)
    return Ellipsoids, Center2, TouchPoint
end




"""
SampleEllipsoids(n, m, p; λ=1.1, γ=1.5)
Creates m ellipsoids in ℝⁿ that intersect. The intersection is regulated by λ.

"""
function SampleEllipsoids(n::Int,  # dimension
    m::Int,  # number of ellipsoids
    p::Real; # sparsity of matrix A
    λ::Real=1.1, # parameter for touching ellipsoid
    γ::Real=1.5 # parameter for making A positive definite
)
    Ellipsoids, CenterEll2, TouchPoint = SampleTwoEllipsoids(n, p, λ=λ, γ=γ)
    point_inter = 0.5 * ((1 - λ)CenterEll2 + (1 + λ)TouchPoint)
    for _ in 3:m
        center = randn(n)
        while any(Ref(center) .∈ Ellipsoids)
            center *= 1.5
        end
        d = λ * (point_inter - center)
        push!(Ellipsoids, GenerateEllipsoid(center, d))
    end
    return Ellipsoids, point_inter
end

"""
    GenerateTouchingEllipsoid(ell, n; λ=1.0)
    Given an ellipsoid ell, generates a touching ellipsoid in ℝⁿ

    """


function GenerateTouchingEllipsoid(ell::EllipsoidCRM,
    n::Int;
    λ::Real=1.1)
    c = randn(n)
    while c ∈ ell
        c *= 1.5
    end
    x̂ = Proj_Ellipsoid(c, ell)
    d = λ * (x̂ - c)
    return GenerateEllipsoid(c, d), c, x̂
end



"""
 GenerateEllipsoid(center, semi_axis)
 Given center and larger semi_axis, generates an ellipsoid in ℝⁿ
"""

function GenerateEllipsoid(center::Vector,
    semi_axis::Vector)
    n = length(center)
    semi_axis_norm = norm(semi_axis)
    Λ = Diagonal([semi_axis_norm; semi_axis_norm .+ 2 * rand(n - 1)])
    # Λ = Diagonal([semi_axis_norm; rand(n - 1)* .8 * semi_axis_norm])   
    Q, _ = qr(randn(n, n))
    M2 = Q' * Λ .^ 2 * Q
    return EllipsoidCRM(center, 0.5(M2 + M2'))
end



##
##################################################################
### Methods
##################################################################


"""
    SucCentCRM_Cyclic(x₀, Ellipsoids)

Successive Centralized Cirumcentered-Reflection method using Cyclic control sequence
"""
function SucCentCRM_Cyclic(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM}; kargs...)
    Projections = Function[x -> Proj_Ellipsoid(x, ell) for ell in Ellipsoids]
    return SucCentCRM_Cyclic(x₀, Projections; kargs...)
end

"""
    SucCentCRM_AlmostViolatedDist(x₀, Ellipsoids)
"""
function SucCentCRM_AlmostViolatedDist(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM}; kargs...)
    Projections = Function[x -> Proj_Ellipsoid(x, ell) for ell in Ellipsoids]
    return SucCentCRM_AlmostViolatedDist(x₀, Projections; kargs...)
end

"""
    SucCentCRM_AlmostViolatedFunc(x₀, Ellipsoids)
"""
function SucCentCRM_AlmostViolatedFunc(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM}; kargs...)
    Projections = Function[x -> Proj_Ellipsoid(x, ell) for ell in Ellipsoids]
    FunctionEval = Function[x -> func_EllipsoidCRM(x, ell) for ell in Ellipsoids]
    return SucCentCRM_AlmostViolatedFunc(x₀, Projections, FunctionEval; kargs...)
end


"""
    SePM(x₀, Ellipsoids)
"""
function SePM(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM}; kargs...)
    Projections = Function[x -> Proj_Ellipsoid(x, ell) for ell in Ellipsoids]
    return SePM(x₀, Projections; kargs...)
end

"""
    CRMProd(x₀, Ellipsoids)
"""
function CRMprod(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM}; kargs...)
    m = length(Ellipsoids)
    n = length(x₀)
    Projections = Function[x -> Proj_Ellipsoid(x, ell) for ell in Ellipsoids]
    Results = CRMprod(x₀, Projections; kargs...)
    Results.iter_total = div(Results.iter_total, 2)
    Results.proj_total = m * n * Results.iter_total
    return Results
end



##
"""
TestEllipsoids()

"""
function TestEllipsoidsRn(
    n::Int, # dimension
    m::Int; # number of ellipsoids
    samples::Int=1,
    λ::Real=1.15,
    ε::Real=1e-6,
    itmax::Int=200_000,
    restarts::Int=1,
    print_file::Bool=false,
    Methods::Vector{Symbol}=[:SucCentCRM_Cyclic, :SucCentCRM_AlmostViolatedDist, :SucCentCRM_AlmostViolatedFunc, :SePM, :CRMprod],
    bench_time::Bool=false,
    gap_distance=false,
    verbose::Bool=false)
    # X = R^n
    # Defines DataFrame for Results
    dfResults, dfFilenames = createDataframes(Methods, projections=true)
    # Fix Random
    Random.seed!(123)
    # Sparsity of first ellipsoid
    p = 2 * inv(n)
    for j = 1:samples
        ℰ, _ = SampleEllipsoids(n, m, p, λ=λ)
        for i = 1:restarts
            x₀ = InitalPoint_EllipsoidCRM(ℰ, n)
            prob_name = savename((Prob=j, Rest=i, dim=n, numsets=m, lambda=λ); equals="", sort=false)
            @info prob_name
            timenow = Dates.now()
            dfrow = []
            dfrowFilename = []
            push!(dfrow, prob_name)
            push!(dfrowFilename, prob_name)
            for mtd in Methods
                func = eval(mtd)
                @info "method: $(string(func))"
                filename = savename("BBILS23", (time=timenow, mtd=mtd, prob=prob_name,), "csv", sort=false)
                print_file ? filedir = datadir("sims", filename) : filedir = ""
                results = func(x₀, ℰ, EPSVAL=ε, filedir=filedir, itmax=itmax, verbose=verbose)
                verbose && @info mtd, results.iter_total, results.proj_total, results.final_tol
                elapsed_time = 0.0
                if bench_time
                    t = @benchmark $func($x₀, $ℰ, EPSVAL=$ε, filedir=$filedir, itmax=$itmax)
                    elapsed_time = (mean(t).time) * 1e-9
                end
                push!(dfrow, results.iter_total)
                push!(dfrow, results.proj_total)
                push!(dfrow, results.final_tol)
                push!(dfrow, elapsed_time)
                push!(dfrowFilename, filedir)
            end
            push!(dfResults, dfrow)
            push!(dfFilenames, dfrowFilename)
        end
    end
    return dfResults, dfFilenames
end
##
#################################################################
### Run tests of Paper
##################################################################
samples = 30
dimensions = [
    (10, 3)
    (10, 10)
    (20, 3)
    (20, 5)
    (20, 10)
    (20, 20)
    (50, 5)
    (50, 10)
    (50, 20)
]
ε = 1e-6
λ = 1.15
itmax = 3_000
Methods = [:SucCentCRM_Cyclic, :SucCentCRM_AlmostViolatedDist, :SucCentCRM_AlmostViolatedFunc, :SePM, :CRMprod]
dfResultsFinal, _ = createDataframes(Methods,projections=true)
for (n, m) ∈ dimensions
    dfResults, _ = TestEllipsoidsRn(n, m, Methods=Methods, bench_time=true, verbose=false, itmax=itmax, samples=samples, λ=λ)
    append!(dfResultsFinal, dfResults)
end

##
#################################################################
### Write results
##################################################################
timenow = Dates.now()
file_name = savename("BBILS23_EllipsoidsCFP", (@dict timenow), "csv")
CSV.write(datadir("sims",file_name), dfResultsFinal)


##
#################################################################
### Make Performance Profiles.
##################################################################
 dfResultsPerprof = CSV.read(datadir("sims",file_name), DataFrame)

 