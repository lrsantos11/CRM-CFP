using LazySets

struct EuclideanBall
    center::Vector{Float64}          # The center of the Ball
    radius::Float64                  # The radius of the ball
    ball::LazySets.Ball2  # LazySets.Ball2 contructor.
    ball_ind::ProximableFunction      # BallL2 translated contructor.
    EuclideanBall(center, radius) = new(center, radius, Ball2(center, radius),
        Translate(IndBallL2(radius), -center))
end

##

function FigureTwoBalls(B1::EuclideanBall,
    B2::EuclideanBall;
    z₀::Vector = [],
    framestyle = :none)
    # Initial Plot 
    plt = plot(size = (400, 400), fg_colour_azis = :lightgray, framestyle = framestyle,
        aspect_ratio = :equal, draw_arrow = true, ticks = :none,
        grid = :none, legend = :topright)
    plot!(B1.ball, c = :dodgerblue)
    plot!(B2.ball, c = :lightgreen)
    # Projection and Reflection  functions
    ProjectB1(z) = ProjectIndicator(B1.ball_ind, z)
    ProjectB2(z) = ProjectIndicator(B2.ball_ind, z)
    ReflectB1(z) = ReflectIndicator(B1.ball_ind, z)
    ReflectB2(z) = ReflectIndicator(B2.ball_ind, z)
    # Plot Initial point
    labelsize = 12
    isempty(z₀) ? z₀ = zeros(2) : nothing
    scatter!(Singleton(z₀), marker = (3, :circle))
    annotate!([(z₀[1] - 0.15, z₀[2] + 0.15, text(L"z", labelsize))])
    zRB1 = ReflectB1(z₀)
    zRB2 = ReflectB2(z₀)
    plot!([Singleton(v) for v in [zRB1, zRB2]], marker = (2, :circle))
    zpCRM = parallelCRMiteration!(z₀, ReflectB1, ReflectB2)
    @show zpCRM ∈ B1.ball
    scatter!(Singleton(zpCRM), marker = (3, :circle), c = :darkorange2)
    annotate!([(zpCRM[1], zpCRM[2] - 0.15, text(L"z_{\mathrm{pCRM}}", labelsize))])
   method_path!(plt, [z₀'; zRB1'], arrow = :none, color = :black, lalpha = 0.2, ls = :dot)
   method_path!(plt, [z₀'; zRB2'], arrow = :none, color = :black, lalpha = 0.2, ls = :dot)
   method_path!(plt, [z₀'; zpCRM'], color = :black, lalpha = 0.5)
    plot!(Ball2(zpCRM, norm(zpCRM - z₀)), alpha = 0.0, ls = :dot, lalpha = 0.2)
    zCent = centralization!(z₀, ProjectB1, ProjectB2)
    @show zCent ∈ B1.ball
    @show zCent ∈ B2.ball
    scatter!(Singleton(zCent), marker = (5, :diamond), c = :bisque4)
    annotate!([(zCent[1] + 0.2, zCent[2], text(L"z_{\mathrm{C}}", labelsize))])
   method_path!(plt, [z₀'; zCent'], color = :black, lalpha = 0.5)
    zC_CRM = parallelCRMiteration!(zCent, ReflectB1, ReflectB2)
    scatter!(Singleton(zC_CRM), marker = (3, :square), c = :red)
    annotate!([(zC_CRM[1] - 0.15, zC_CRM[2] - 0.15, text(L"z_{\mathrm{cCRM}}", :red, labelsize))])
   method_path!(plt, [zCent'; zC_CRM'], color = :black, lalpha = 0.5)
    @show zC_CRM ∈ B1.ball
    @show zC_CRM ∈ B2.ball
    return plt, zpCRM, zCent, zC_CRM
end
##
labelsize = 12
B1 = EuclideanBall([0, 1.0], 2.0)
B2 = EuclideanBall([0, -5.0], 4.0)
z₀ = [2.5, 2.2]
framestyle = :none
plt1, zpCRM, zCent, zC_CRM = FigureTwoBalls(B1, B2, z₀ = z₀, framestyle = framestyle)
xlims!(plt1, -0.5, 4.5, aspect_ratio = :equal)
ylims!(plt1, -2.5, 2.5, aspect_ratio = :equal)
annotate!([(0.0, 1.5, text(L"X", labelsize))])
annotate!([(0, -2, text(L"Y", labelsize))])

##

savefig(plt1, plotsdir("CentralizedvsNonCentralized.pdf"))

##
# z₀ = [6.2, 0.25]
# plt2,  = FigureTwoBalls(B1,B2,z₀ = z₀, framestyle = :none)
# xlims!(plt2, -3.5,7.5,aspect_ratio = :equal)
# ylims!(plt2, -7.5,3.5,aspect_ratio = :equal)
labelsize = 12
plt2 = plot(framestyle = :none, aspect_ratio = :equal, size = (400, 400))
plot!(B1.ball, c = :dodgerblue)
plot!(B2.ball, c = :lightgreen)
ProjectB1(z) = ProjectIndicator(B1.ball_ind, z)
ProjectB2(z) = ProjectIndicator(B2.ball_ind, z)
ReflectB1(z) = ReflectIndicator(B1.ball_ind, z)
ReflectB2(z) = ReflectIndicator(B2.ball_ind, z)
zPB1 = ProjectB1(zCent)
zPB2 = ProjectB2(zCent)
zRB1 = ReflectB1(zCent)
zRB2 = ReflectB2(zCent)
plot!([Singleton(v) for v in [zPB1, zPB2, zRB1, zRB2]], marker = (2, :circle))
xlims!(plt2, 0, 1.0, aspect_ratio = :equal)
ylims!(plt2, -1.5, -0.5, aspect_ratio = :equal)
scatter!(Singleton(zCent), marker = (5, :diamond), c = :bisque4)
annotate!([(zCent[1] - 0.05, zCent[2] - 0.05, text(L"z_{\mathrm{C}}", labelsize))])
scatter!(Singleton(zC_CRM), marker = (3, :square), c = :red)
annotate!([(zC_CRM[1], zC_CRM[2] - 0.05, text(L"z_{\mathrm{cCRM}}", :red, labelsize))])
plot!(Ball2(zC_CRM, norm(zC_CRM - zCent)), alpha = 0.0, ls = :dot, lalpha = 0.2)
MethodPath!(plt2, [zCent'; zC_CRM'], color = :black, lalpha = 0.5, λ = 0.88)
MethodPath!(plt2, [z₀';  zCent'], color = :black, lalpha = 0.5, λ = 0.99)
MethodPath!(plt2, [zCent'; zRB1'], arrow = :none, color = :black, lalpha = 0.2, ls = :dot)
MethodPath!(plt2, [zCent'; zRB2'], arrow = :none, color = :black, lalpha = 0.2, ls = :dot)
v₁ = zCent - zPB1
b₁ = dot(v₁, zPB1)
plot!(x -> (b₁ - v₁[1] * x) / v₁[2], -0.2, 1, leg = false, color = :dodgerblue, lalpha = 0.8, ls = :dashdot)
v₂ = zCent - zPB2
b₂ = dot(v₂, zPB2)
plot!(x -> (b₂ - v₂[1] * x) / v₂[2], -0.2, 1, leg = false, color = :green, lalpha = 0.8, ls = :dashdot)
annotate!([(0.3, -0.7, text(L"X", labelsize))])
annotate!([(0.3, -1.3, text(L"Y", labelsize))])
annotate!([(0.15, -1.1, text(L"S_X", :dodgerblue, labelsize))])
annotate!([(0.15, -.96, text(L"S_Y", :green, labelsize))])

##

savefig(plt2, plotsdir("Centralized_Halfspaces.pdf"))
