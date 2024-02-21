using LazySets

struct EuclideanBall
    center::Vector{Float64}          # The center of the Ball
    radius::Float64                  # The radius of the ball
    ball::LazySets.Ball2  # LazySets.Ball2 contructor.
    ball_ind::Translate     # BallL2 translated contructor.
    EuclideanBall(center, radius) = new(center, radius, Ball2(center, radius),
        Translate(IndBallL2(radius), -center))
end

function PACA(x₀::AbstractVector,
    Balls::AbstractVector{EuclideanBall},
    ϵ::Function;
    kargs...)
    FuncBalls = Function[x -> g_ball(x, ball) for ball in Balls]
    ∂Balls = Function[x -> ∂g_ball(x, ball) for ball in Balls]
    return PACA(x₀, FuncBalls, ∂Balls, ϵ; kargs...)
end

ϵ(k) = 1 / (sqrt(k + 1))

g_ball(X, ball) = (X[1] - ball.center[1])^2 + (X[2] - ball.center[1]) - ball.radius
∂g_ball(X, ball) = [2 * (X[1] - ball.center[1]), 2 * (X[2] - ball.center[2])]

##

# framestyle = :origin
framestyle = :none
labelsize = 16
B1 = EuclideanBall([0, 0.0], 4.5)
B2 = EuclideanBall([1.5, 0.], 1.6)
B3 = EuclideanBall([1, 0.], 1.3)
B4 = EuclideanBall([.5, 0.0], 1.1)
Balls = [B1, B2, B3, B4]

z₀ = [5.5, 1.0]
PACAfile = datadir("sims", "PACA_TwoBalls.csv")

resultsPACA = PACA(z₀, Balls, ϵ, verbose=true, filedir=PACAfile)

xPACA = readdlm(PACAfile)
PACA_iter_total = Int(xPACA[end, 1])
framestyle = :none


##
# function FigureTwoBalls(B1::EuclideanBall,
#     B2::EuclideanBall;
#     z₀::Vector = [],
#     framestyle = :none)
#     # Initial Plot 
    plt = plot(size = (900, 760), fg_colour_axis = :lightgray, framestyle = framestyle,
        aspect_ratio = :equal, draw_arrow = true, 
        ticks = :none,
        grid = :none, legend = :topright)
    plot!(B1.ball, c = :dodgerblue)
    annotate!([(3.4, -0.5, text(L"C", labelsize))])
    plot!(B2.ball, c = :lightgreen)
    annotate!([(2.4, -0.5, text(L"C^{k}", labelsize))])
    plot!(B3.ball, c=:yellow)
    annotate!([(1.8, -0.5, text(L"C^{k-1}", labelsize))])
    plot!(B4.ball, c=:darkgreen)
    annotate!([(1.235, -0.5, text(L"C^{1}", labelsize))])
annotate!([(1.45, -0.52, text(L"\cdots", labelsize))])

xlims!(plt, 1, 6, aspect_ratio=:equal)

ylims!(plt, -2.0, 2.0, aspect_ratio=:equal)





    scatter!(Singleton(z₀), marker=(4, :circle))
points = [
            [5.3,0.9],
            [4.6, .45 ],
            [4.1,0.4],

]
label_points!(plt, [z₀], label_size=labelsize, xshift=0.12, yshift=-0.01, var_name="x^{")
# MethodPath!(plt, xPACA[1:5, 3:4], c=:dodgerblue2, alpha=0.4, λ=0.95)

    plot!(plt, [Singleton(v) for v in points], c=:dodgerblue2, alpha=0.8, ms=6, m=:diamond)
label_points!(plt, [points[1]], label_size=labelsize, xshift=0.12, yshift=-0.1, var_name="x^{", var_indices=[1])
label_points!(plt, [points[2]], label_size=labelsize, xshift=0.15, yshift=-0.1, var_name="x^{k-", var_indices=[1])
annotate!([(points[3][1] + 0.12, points[3][2] - 0.1, text(L"x^k", labelsize))])

MethodPath!(plt, Matrix([z₀ points[1]]'), c=:dodgerblue2, alpha=0.4, λ=0.95)

MethodPath!(plt, Matrix([points[2] points[3]]'), c=:dodgerblue2, alpha=0.4, λ=0.95)



# label_points!(plt, [v for v in points], num_points=3, label_size=11, xshift=0.1, yshift=-0.15, var_name="x^{k-", var_indices = [1,2,3])
# scatter!(Singleton(resultsPACA.xApprox), marker=(4, :circle), c=:red)

# annotate!([(resultsPACA.xApprox[1] - 0.2, resultsPACA.xApprox[2] -.12 , text(L"\bar{z}", 14))])

##


##
savefig(plt, plotsdir("PACA_Figure.pdf"))
# cp(plotsdir("PACA_Figure1.pdf"), "../paper/figures/PACA_Figure1.pdf")
plt

