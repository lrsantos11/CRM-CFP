
##################################################################
## Creating Samples
##################################################################


function createTwoEllipsoids(n::Int,
    p::Real;
    λ::Real = 1.0)
    Ellipsoids = EllipsoidCRM[]
    A = Matrix(sprandn(n, n, p))
    γ = 1.5
    A = (γ * I + A' * A)
    a = rand(n)
    b = A * a
    adotAa = dot(a, b)
    b .*= -1.0
    α = (1 + γ) * adotAa
    push!(Ellipsoids, EllipsoidCRM(A, b, α))
    TouchEll, TouchPoint = TouchingEllipsoid(Ellipsoids[1], n, λ = λ)
    push!(Ellipsoids, TouchEll)
    return Ellipsoids, TouchPoint
end

function TouchingEllipsoid(ell::EllipsoidCRM,
    n::Int;
    λ::Real = 1.0)
    c = randn(n)
    while c ∈ ell
        c *= 1.5
    end
    x̂ = Proj_Ellipsoid(c, ell)
    d = λ * (x̂ - c)
    Λ = Diagonal([norm(d); norm(d) .+ 2 * rand(n - 1)])
    Q, = qr(d)
    A2 = Λ * Q
    M2 = Symmetric(A2' * A2)
    return EllipsoidCRM(c, M2), x̂
end


function createEllipsoids(n::Int, p::Real, m::Int)
    Ellipsoids = EllipsoidCRM[]
    for index = 1:m
        A = sprandn(n, n, p)
        γ = 1.5
        A = (γ * I + A' * A)
        a = rand(n)
        b = A * a
        adotAa = dot(a, b)
        b .*= -1.0
        α = (1 + γ) * adotAa
        push!(Ellipsoids, EllipsoidCRM(A, b, α))
    end
    return Ellipsoids
end

##################################################################
## Methods
##################################################################

"""
centralizedCRM(x₀, Ellipsoids)
Uses cCRM to find a point into intersection of two EllipsoidCRM 
"""
function centralizedCRM(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM};
    kwargs...)
    ProjectA(x) = Proj_Ellipsoid(x, Ellipsoids[1])
    ProjectB(x) = Proj_Ellipsoid(x, Ellipsoids[2])
    return centralizedCRM(x₀, ProjectA, ProjectB; kwargs...)
end


"""
CRMprod(x₀, Ellipsoids)
Uses cCRM to find a point into intersection of two EllipsoidCRM 
"""
function CRMprod(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM};
    kwargs...)
    ProjectA(x) = Proj_Ellipsoid(x, Ellipsoids[1])
    ProjectB(x) = Proj_Ellipsoid(x, Ellipsoids[2])
    return CRMprod(x₀, ProjectA, ProjectB; kwargs...)
end



"""
MAP(x₀, Ellipsoids)
Uses MAP to find  a point into intersection of two EllipsoidCRM 
"""
function MAP(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM};
    kwargs...)
    ProjectA(x) = Proj_Ellipsoid(x, Ellipsoids[1])
    ProjectB(x) = Proj_Ellipsoid(x, Ellipsoids[2])
    return MAP(x₀, ProjectA, ProjectB; kwargs...)
end


"""
SPM(x₀, Ellipsoids)
Uses SPM to find  a point into intersection of two EllipsoidCRM 
"""
function SPM(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM};
    kwargs...)
    ProjectA(x) = Proj_Ellipsoid(x, Ellipsoids[1])
    ProjectB(x) = Proj_Ellipsoid(x, Ellipsoids[2])
    return SPM(x₀, ProjectA, ProjectB; kwargs...)
end


"""
DRM(x₀, Ellipsoids)
Uses DRM to find  a point into intersection of two EllipsoidCRM 
"""
function DRM(x₀::Vector,
    Ellipsoids::Vector{EllipsoidCRM};
    kwargs...)
    ProjectA(x) = Proj_Ellipsoid(x, Ellipsoids[1])
    ProjectB(x) = Proj_Ellipsoid(x, Ellipsoids[2])
    return DRM(x₀, ProjectA, ProjectB; kwargs...)
end


##################################################################
## Generating 2D plots
##################################################################
"""
makeplots_Ellipsoids2D()
Makes plots for Ellipsoids in 2D
"""
function makeplots_TwoEllipsoids2D(λ::Real;
    ε::Real=1e-4,
    itmax::Int=10_000,
    generate_results::Bool=false)
    ##
    Random.seed!(11)
    n = 2
    ℰ, xSol = createTwoEllipsoids(n, 2 / n, λ=λ)
    x₀ = [0.5, 2.8]
    while x₀ ∈ ℰ[1] || x₀ ∈ ℰ[2]
        x₀ *= 1.2
    end
    plt = plot(Ellipsoid.(ℰ),
        framestyle=:none,
        leg=:bottomright, fillcolor=[:turquoise4 :palegreen4],
        aspect_ratio=:equal,
        legendfontsize=14
    )
    MAPfile = datadir("sims", savename("BBIS21_TwoEllip2DMAP", (@dict λ), "csv"))
    SPMfile = datadir("sims", savename("BBIS21_TwoEllip2DSPM", (@dict λ), "csv"))
    cCRMfile = datadir("sims", savename("BBIS21_TwoEllip2DcCRM", (@dict λ), "csv"))
    CRMprodfile = datadir("sims", savename("BBIS21_TwoEllip2DCRMprod", (@dict λ), "csv"))
    if generate_results
        resultsMAP = MAP(x₀, ℰ, filedir=MAPfile, itmax=itmax, EPSVAL=ε)
        resultsSPM = SPM(x₀, ℰ, filedir=SPMfile, itmax=itmax, EPSVAL=ε)
        resultscCRM = centralizedCRM(x₀, ℰ, filedir=cCRMfile, itmax=itmax, EPSVAL=ε)
        resultsCRMprod = CRMprod(x₀, ℰ, filedir=CRMprodfile, itmax=itmax, EPSVAL=ε)
    end
    ##
    xMAP = readdlm(MAPfile)
    xSPM = readdlm(SPMfile)
    xcCRM = readdlm(cCRMfile)
    xCRMprod = readdlm(CRMprodfile)
    MAP_iter_total = Int(xMAP[end, 1])
    SPM_iter_total = Int(xSPM[end, 1])
    cCRM_iter_total = Int(xcCRM[end, 1])
    CRMprod_iter_total = Int(xCRMprod[end, 1])

    ##
    plot!([Singleton(v) for v in eachrow(xSPM[2:end, 3:4])], label="SPM -- $(SPM_iter_total) projections", c=:dodgerblue2, alpha=0.8, ms=4, m=:utriangle)
    plot!([Singleton(v) for v in eachrow(xCRMprod[1:end, 3:4])], label="CRMprod -- $(CRMprod_iter_total) projections", c=:yellow, alpha=0.8, ms=3, m=:square)
    plot!([Singleton(v) for v in eachrow(xMAP[3:2:end, 3:4])], label="MAP -- $(MAP_iter_total) projections", c=:red, alpha=0.8, ms=3)
    # plot!([Singleton(v) for v in eachrow(xMAP[2:2:100,3:4])],c = :tomato, alpha = 0.4, ms = 3)
    plot!([Singleton(v) for v in eachrow(xcCRM[2:end, 3:4])], label="cCRM -- $(cCRM_iter_total) projections", c=:green, alpha=0.8, ms=4, m=:diamond)
    MethodPath!(plt, xMAP[1:2:9, 3:4], c=:red)
    MethodPath!(plt, xSPM[1:8, 3:4], c=:dodgerblue2)
    MethodPath!(plt, xcCRM[:, 3:4], c=:green)
    MethodPath!(plt, xCRMprod[1:4, 3:4], c = :yellow)
    plot!(Singleton(x₀), c=:black, m=:square, alpha=1)
    label_points!(plt, [x₀], label_size=14, xshift=0.25, yshift=0.01)
    return plt
end
##

λ = 1.14
plt1 = makeplots_TwoEllipsoids2D(λ, itmax = 1_000, generate_results = true)
figname = savename("BBIS21_Ellipsoids", (@dict λ), "pdf")
savefig(plt1, plotsdir(figname))
plt1

##

λ = 1.0
plt2 = makeplots_TwoEllipsoids2D(λ, itmax = 100_000, generate_results = false)
figname = savename("BBIS21_Ellipsoids", (@dict λ), "pdf")
savefig(plt2, plotsdir(figname))
plt2
##################################################################
## Generating Performance Profiles 
###################################################################

"""
TestEllipsoids()

"""
function tests_TwoEllipsoidsRn(; n::Int = 100,
    samples::Int = 1,
    λ::Real = 1.0,
    ε::Real = 1e-6,
    itmax::Int = 200_000,
    restarts::Int = 1,
    print_file::Bool = false,
    Methods::Vector{Symbol} = [:centralizedCRM, :DRM],
    bench_time::Bool = false,
    gap_distance = false)
    # X = R^n
    # Defines DataFrame for Results
    dfResults, dfFilenames = createDataframes(Methods)
    # Fix Random
    Random.seed!(1234)
    # Sparsity of first ellipsoid
    p = 2 * inv(n)
    for j = 1:samples
        ℰ, xSol = createTwoEllipsoids(n, p, λ = λ)
        for i = 1:restarts
            x₀ = StartingPoint(n)
            while x₀ ∈ ℰ[1] || x₀ ∈ ℰ[2]
                x₀ *= 1.2
            end
            prob_name = savename((Prob = j, Rest = i, n = n, λ = λ))
            @show prob_name
            timenow = Dates.now()
            dfrow = []
            dfrowFilename = []
            push!(dfrow, prob_name)
            push!(dfrowFilename, prob_name)
            for mtd in Methods
                func = eval(mtd)
                filename = savename("BBIS21", (mtd = mtd, time = timenow, prob = prob_name), "csv", sort = false)
                print_file ? filedir = datadir("sims", filename) : filedir = ""
                results = func(x₀, ℰ, EPSVAL = ε, gap_distance = gap_distance, filedir = filedir, itmax = itmax, xSol = xSol)
                @show mtd, results.iter_total, results.final_tol
                elapsed_time = 0.0
                if bench_time
                    t = @benchmark $func($x₀, $ℰ, $itmax, EPSVAL = $ε, gap_distance = gap_distance, filedir = $filedir, xSol = $xSol)
                    elapsed_time = (mean(t).time) * 1e-9
                end
                push!(dfrow, results.iter_total)
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

# Slater point in the intersection
ε = 1e-6
samples = 30
Methods = [:CRMprod]
λ = 1.1
itmax = 10_000
# dfResul/ts, dfFilesname = tests_TwoEllipsoidsRn(samples = samples, ε = ε, Methods = Methods, λ = λ, itmax = itmax, gap_distance = true)
##
## To write data. 
# timenow = Dates.now()
# CSV.write(datadir("sims", savename("BBIS21_TwoEllipsoidsTableSlaterPoint", (@dict timenow λ), "csv")), dfResultsTotal)

## # To make Performance Profiles.
 dfResultsTotal = CSV.read(datadir("sims","BBIS21_TwoEllipsoidsTableSlaterPoint_timenow=2021-11-09T10:15:08.589_λ=1.1.csv"), DataFrame)

T = Float64[dfResultsTotal.centralizedCRM_it dfResultsTotal.MAP_it dfResultsTotal.CRMprod_it]
T[findall(row -> row >= 10000, T)] .= Inf
perprof1 = performance_profile(PlotsBackend(),
    T, 
    # logscale = false,
    ["cCRM", "MAP", "CRMprod"],
    legend=:bottomright, framestyle=:box,
    linestyles=[:solid, :dash, :dot])
ylabel!("Fraction of problems solved")
xticks!(perprof1,[0, 2, 4, 6], [L"2^0", L"2^2", L"2^4", L"2^6"])
# title!(L"Performance Profile -- Total projections comparison -- tolerance $\varepsilon = 10^{-4}$")
savefig(perprof1, plotsdir(savename("BBIS21_TwoEllipsoids_Perprof", (@dict λ), "pdf")))
perprof1
##
# Summa = describe(dfResults[:, [:centralizedCRM_it, :MAP_it]], :mean, :std, :median, :min, :max)


## No Slater point in the intersectio = [:centralizedCRM, :MAP]
λ = 1.0
samples = 15
dfResults, dfFilesname = tests_TwoEllipsoidsRn(samples = samples, ε = 1e-6, Methods = Methods, λ = λ, itmax = 500_000, gap_distance = false)
## To write data. 
timenow = Dates.now()
CSV.write(datadir("sims", savename("BBIS21_TwoEllipsoidsTableNoSlaterPoint", (@dict timenow λ), "csv")), dfResults)



##
perprof2 = performance_profile(PlotsBackend(),
    Float64.(hcat(dfResults.centralizedCRM_it, dfResults.MAP_it)),
    ["cCRM", "MAP"],
    legend = :bottomright,
    framestyle = :box,
    linestyles = [:solid, :dash])
ylabel!("Percentage of problems solved")
# title!(L"Performance Profile -- Total projections comparison -- tolerance $\varepsilon = 10^{-4}$")
savefig(perprof2, plotsdir(savename("BBIS21_TwoEllipsoids_Perprof", (@dict λ), "pdf")))
perprof2


## To consult the results.

# dfResultsEllips = CSV.read(datadir("sims","BBIS21_TwoEllipsoidsTable.csv"), DataFrame)
# @show df_Summary = describe(dfResultsEllips,:mean,:max,:min,:std,:median)[2:end,:]

##
MersenneTwister(123)
n = 2
ℰ, xSol = createTwoEllipsoids(n, 2 / n, λ = λ)
x₀ = [0.5, 2.8]
while x₀ ∈ ℰ[1] || x₀ ∈ ℰ[2]
    x₀ *= 1.2
end
plt = plot(Ellipsoid.(ℰ),
    framestyle = :none,
    leg = :bottomleft, fillcolor = [:turquoise4 :palegreen4],
    aspect_ratio = :equal,
    legendfontsize = 14
)