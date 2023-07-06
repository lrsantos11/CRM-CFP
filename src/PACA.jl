"""
    PACA(x₀, Projections)

Perturbed Approximate Circumcenter Algorithm
"""
function PACA(x₀::Vector,
              Functions::Vector{Function},
              Subgrads::Vector{Function},
              ϵ::Function;
                EPSVAL::Float64=1e-12, 
                itmax::Int=100,
                filedir::String="",
                verbose::Bool=false,
)
    xPACA = copy(x₀)
    wₖ = similar(x₀)
    m = length(Functions)
    invm = inv(m)
    vₖ = [copy(x₀) for _ in 1:m]
    k = 0
    tol = 1.0
    printOnFile(filedir, k, tol, xPACA, deletefile=true)
    solved = false
    tired = false
    while !(solved || tired)
        ϵₖ = ϵ(k)
        computevₖ!(vₖ, xPACA, Functions, Subgrads, m, ϵ=ϵₖ)
        copyto!(wₖ, mean(vₖ))
        αₖ = (mapreduce(v -> dot(v, v), +, vₖ) * invm) / dot(wₖ, wₖ)
        wₖ *= αₖ
        xPACA .-= wₖ
        k += 1
        tol = maximum([f(xPACA) for f in Functions])
        verbose && @info "tol = $tol"
        solved = tol ≤ EPSVAL
        tired = k ≥ itmax
        printOnFile(filedir, k, tol, xPACA)
    end
    return Results(iter_total=k, proj_total= m*k, final_tol=tol, xApprox=xPACA, method=:PACA)
end

    
"""
    MCSPM(x₀, Projections)

Modified Cyclic Subgradient Projection method  from [Pierro:1988]

[Pierro:1986] A. R. De Pierro and A. N. Iusem, “A finitely convergent ‘row-action’ method for the convex feasibility problem,” Appl Math Optim, vol. 17, no. 1, Art. no. 1, Jan. 1988, doi: 10.1007/BF01448368.

"""
function MCSPM(x₀::Vector,
    Functions::Vector{Function},
    Subgrads::Vector{Function},
    ϵ::Function;
    EPSVAL::Float64=1e-12,
    itmax::Int=100,
    filedir::String="",
    verbose::Bool=false
)
    xMCSPM = copy(x₀)
    m = length(Functions)
    @assert m == length(Subgrads)
    k = 0
    tol = 1.0
    printOnFile(filedir, k, tol, xMCSPM, deletefile=true)
    solved = false
    tired = false
    while !(solved || tired)
        ϵₖ = ϵ(k)
        for i in 1:m
            @inbounds xMCSPM .-= computevₖⁱ(xMCSPM, Functions[i], Subgrads[i], ϵ=ϵₖ)
        end
        k += 1
        tol = maximum([f(xMCSPM) for f in Functions])
        verbose && @info "tol = $tol"
        solved = tol ≤ EPSVAL
        tired = k ≥ itmax
        printOnFile(filedir, k, tol, xMCSPM)
    end
    return Results(iter_total=k, proj_total=m * k, final_tol=tol, xApprox=xMCSPM, method=:MCSPM)
end


"""
    MSSPM(x₀, Projections, Subgradients)
Modified simultaneous Subgradient Projection method for ellipsoids from [Iusem:1986]

[Iusem:1986] A. N. Iusem and L. Moledo, “A finitely convergent method of simultaneous subgradient projections for the convex feasibility problem,” Matemática Aplicada e Computacional, vol. 5, no. 2, pp. 169–184, 1986.

"""
function MSSPM(x₀::Vector,
    Functions::Vector{Function},
    Subgrads::Vector{Function},
    ϵ::Function;
    EPSVAL::Float64=1e-12,
    itmax::Int=100,
    filedir::String="",
    verbose::Bool=false
)
    xMSSPM = copy(x₀)
    m = length(Functions)
    @assert m == length(Subgrads)
    vₖ = [copy(x₀) for _ in 1:m]
    k = 0
    tol = 1.0
    printOnFile(filedir, k, tol, xMSSPM, deletefile=true)
    solved = false
    tired = false
    while !(solved || tired)
        ϵₖ = ϵ(k)
        computevₖ!(vₖ, xMSSPM, Functions, Subgrads, m, ϵ=ϵₖ)
        xMSSPM .-= mean(vₖ)
        k += 1
        tol = maximum([f(xMSSPM) for f in Functions])
        verbose && @info "tol = $tol"
        solved = tol ≤ EPSVAL
        tired = k ≥ itmax
        printOnFile(filedir, k, tol, xMSSPM)
    end
    return Results(iter_total=k, proj_total=m * k, final_tol=tol, xApprox=xMSSPM, method=:MSSPM)
end




"""
    computevₖ(x₀, func_f, ∂f)

"""
function computevₖⁱ(x::AbstractVector, func_f::Function, ∂f::Function;
                  ϵ::Real = 0.0 # perturbation
                  )
    fx = func_f(x)
    ∂fx = ∂f(x)
    return  (max(0.0, fx + ϵ) / dot(∂fx, ∂fx)) * ∂fx
end
    
    """
        computevₖ!(vₖ, x₀, Functions, Subgrads)
    
    """
function computevₖ!(vₖ, x, Functions, Subgrads, m;
                  ϵ::Number = 0.0 # perturbation
                  )
    Threads.@threads for i in 1:m
        copyto!(vₖ[i], computevₖⁱ(x, Functions[i], Subgrads[i], ϵ=ϵ))
    end
end