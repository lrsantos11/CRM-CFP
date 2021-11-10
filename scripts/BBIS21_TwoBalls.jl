##################################################################
## Basic Functions for Two Balls tests
##################################################################
using BenchmarkTools
Random.seed!(1)
n = 100
T = BigFloat
xSol = zeros(T,n)
center1 = T.(StartingPoint(n))
center2 = -100*center1
𝓑₁ = Translate(IndBallL2(norm(center1)), center1)
𝓑₂ = Translate(IndBallL2(norm(center2)), center2)
x₀ = T.(StartingPoint(n))
itmax = 1000000
Proj_BallA(x::AbstractVector) = ProjectIndicator(𝓑₁, x)
Proj_BallB(x::AbstractVector) = ProjectIndicator(𝓑₂, x)
ε = 1e-16
##
centralizedCRM(x₀, Proj_BallA, Proj_BallB, xSol = xSol, gap_distance = false, EPSVAL = ε, itmax = itmax)
##
MAP(x₀, Proj_BallA, Proj_BallB, xSol = xSol, gap_distance = false, EPSVAL = ε, itmax = itmax)
##
DRM(x₀, Proj_BallA, Proj_BallB, xSol = xSol, gap_distance = true, EPSVAL = ε, itmax = itmax)





