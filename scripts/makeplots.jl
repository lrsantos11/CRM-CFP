include("auxplots.jl")
#
itmax = 10^6
ε = 1e-3
x₀ = [1.1,0]
# Quadratic without error bound
α, β = 1.0, 0.0
xSol = [0.0, 0.0]
plt1 = FigureQuadratic(α,β,x₀,xSol,ε,itmax)
savefig(plt1,plotsdir("FigureQuadraticNonerrorbound.pdf"))
##
# Quadratic without error bound
α, β = 1.0, -0.16
xSol = [0.4, 0.0]
plt2 = FigureQuadratic(α,β,x₀,xSol,ε,itmax)
savefig(plt2,plotsdir("FigureQuadraticWitherrorbound.pdf"))
