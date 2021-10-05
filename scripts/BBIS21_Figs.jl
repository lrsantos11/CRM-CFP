using LazySets

include(scriptsdir("plots_util.jl"))
struct EuclideanBall
    center :: Vector{Float64}          # The center of the Ball
    radius :: Float64                  # The radius of the ball
    ball :: LazySets.Ball2  # LazySets.Ball2 contructor.
    ball_ind :: ProximableFunction      # BallL2 translated contructor.
    EuclideanBall(center, radius) = new(center, radius, Ball2(center, radius),
                                         Translate(IndBallL2(radius),-center))
end

##

function FigureTwoBalls(B1::EuclideanBall, 
                        B2::EuclideanBall; 
                        z₀::Vector = [], 
                        framestyle = :none)
    # Initial Plot 
    plt = plot(size=(400,400), fg_colour_azis=:lightgray, framestyle=framestyle,
                aspect_ratio=:equal, draw_arrow=true, ticks=:none,
                grid=:none, legend=:topright)
    plot!(B1.ball, c = :dodgerblue)  
    plot!(B2.ball, c = :lightgreen) 
    # Projection and Reflection  functions
    ProjectB1(z) =  ProjectIndicator(B1.ball_ind,z)  
    ProjectB2(z) =  ProjectIndicator(B2.ball_ind,z)
    ReflectB1(z) =  ReflectIndicator(B1.ball_ind,z)  
    ReflectB2(z) =  ReflectIndicator(B2.ball_ind,z)    
    # Plot Initial point
    isempty(z₀) ? z₀ = zeros(2) : nothing
    scatter!(Singleton(z₀),marker=(3,:circle)) 
    zRB1 = ReflectB1(z₀)
    zRB2 = ReflectB2(z₀)
    plot!([Singleton(v) for v in [zRB1, zRB2]],marker=(2,:circle))
    zpCRM = parallelCRMiteration!(z₀, ReflectB1, ReflectB2)
    @show zpCRM ∈ B1.ball
    scatter!(Singleton(zpCRM),marker=(3,:circle))
    MethodPath!(plt, [z₀' ; zRB1'; zRB2'; z₀' ], arrow=:none)
    MethodPath!(plt, [z₀' ; zpCRM' ], color = :yellow)
    plot!(Ball2(zpCRM,norm(zpCRM-z₀)),alpha = 0.1,ls=:dot,lalpha=0.2)
    zCent = centralization!(z₀, ProjectB1, ProjectB2)
    @show zCent ∈ B1.ball
    @show zCent ∈ B2.ball
    scatter!(Singleton(zCent),marker=(2,:circle))
    zC_CRM = parallelCRMiteration!(zCent, ReflectB1, ReflectB2)
    scatter!(Singleton(zC_CRM),marker=(2,:circle), c = :red)
    MethodPath!(plt, [z₀' ; zCent'; zC_CRM' ], color = :red)
    @show zC_CRM ∈ B1.ball
    @show zC_CRM ∈ B2.ball
    return plt, zpCRM
end
##
B1 =  EuclideanBall([0,1.],2.)
B2 =  EuclideanBall([0,-5.],4.)
z₀ = [2.5, 2.2]
framestyle = :none
plt1, zpCRM = FigureTwoBalls(B1, B2, z₀ = z₀, framestyle = framestyle)
xlims!(plt1, -0.5,4.5,aspect_ratio = :equal)
ylims!(plt1, -2.5,2.5,aspect_ratio = :equal)

##
z₀ = [6.2, 0.25]
plt2,  = FigureTwoBalls(B1,B2,z₀ = z₀, framestyle = :none)
xlims!(plt2, -3.5,7.5,aspect_ratio = :equal)
ylims!(plt2, -7.5,3.5,aspect_ratio = :equal)
##
plot(size=(800,400),plt1,plt2, layout = (1,2))
##

savefig(plotsdir("CentralizedvsNonCentralized.pdf"))

##
# """
#     FigureQuadraticApprox(α,β,x₀,xSol,ε)
# """
# function FigureQuadraticApprox(α::Number,β::Number, x₀::Vector, xSol::Vector,ε::Number; itmax::Int64=6,  
#     methods = [:MAP, :CRM], mrk_color = [:green, :blue], max_iter_plotted ::Int=13, ymax::Number = 1.25)
#     # On the Right
#     # f(t) = α t^2 + β 
#     shiftvec = [0,β]
#     Hyper = IndAffine([0,1.], 0.0)
#     g(x) = α*dot(x[1:end-1],x[1:end-1]) - x[end]
#     ∂g(x) = [2*α*x[1:end-1]; -1]
#     ProjectARight(x) =  ApproxProject(x-shiftvec,g,∂g) + shiftvec
#     ProjectBRight(x) =  ProjectIndicator(Hyper,x)
#     plt = plot(size=(400,400),fg_colour_axis=:lightgray,framestyle=:origin,
#                aspect_ratio=:equal,draw_arrow=true,ticks=:none,grid=:none,legend=:topright)
#     plot!(x -> α*x^2 + β,-0.5,sqrt((ymax - β)/α),lw=1,c=:blue,fill = (ymax, 0.5, :dodgerblue),label="")
#     scatter!([x₀[1]],[x₀[2]],label="",marker=(3,:circle))
#     timenow = Dates.now()
#     filedirs = String[]
#     for (index,mtd) in enumerate(methods)
#         filename = savename("AABBIS21FigApprox_Quadratic",(mtd=mtd,date=timenow),"csv",sort=false)
#         filedir = datadir("sims",filename)
#         push!(filedirs,filedir)
#         func = getfield(Main,mtd)
#         Result = func(x₀, ProjectARight,ProjectBRight,itmax=itmax,filedir=filedir,EPSVAL=ε,xSol=xSol,print_intermediate = true)
#         @show Result
#         pts_read_mtd = 2*min(Result.iter_total,max_iter_plotted) + 1
#         xmtd = (readdlm(filedir))[1:pts_read_mtd,3:end]
#         scatter_points_Proj_x = @view xmtd[2:2:end,1]
#         scatter_points_Proj_y = @view xmtd[2:2:end,2]
#         scatter_points_mtd_x = @view xmtd[3:2:end,1]
#         scatter_points_mtd_y = @view xmtd[3:2:end,2]
#         scatter!(scatter_points_Proj_x,scatter_points_Proj_y, c=:blue, marker=(2,:circle),
#                 label="")                
#         scatter!(scatter_points_mtd_x,scatter_points_mtd_y, c=mrk_color[index], marker=(3,:circle),
#                 label="A$(String(mtd)) ($(Result.iter_total) it.)")        
#         MethodPath(plt,xmtd)
#     end
#     return plt,filedirs
# end


# #
# itmax = 10^6
# ε = 1e-3
# x₀ = [1.1,0]
# ##
# # Quadratic without error bound
# α, β = 1., 0.0
# xSol = [0.0, 0.0]
# # 
# mtd = :CRM
# pltCARM1,filedirs = FigureQuadraticApprox(α,β,x₀,xSol,ε,itmax=itmax,methods=[mtd],mrk_color = [:red], max_iter_plotted=5)
# annotate!([(0.35, .75 ,text(L"\phi(t) = |t|^\alpha",12))])
# label_points!(pltCARM1,filedirs[1],7,print_proj=true,num_print_proj=2,ProjName="P^S")
# ##
# mtd = :MAP
# pltMAAP1,filedirs = FigureQuadraticApprox(α,β,x₀,xSol,ε,itmax=itmax,methods=[mtd], mrk_color = [:green], max_iter_plotted=21)
# annotate!([(0.35, .75 ,text(L"\phi(t) = |t|^\alpha",12))])
# label_points!(pltMAAP1,filedirs[1],7,print_proj=true,var_name="z",num_print_proj=2,ProjName="P^S")

# ##
# # Quadratic with error bound \phi(t) = αx² + β 

# α, β = 1., -0.06
# xSol = [sqrt(abs(β)/α), 0.0]
# ##
# mtd = :CRM
# pltCARM2,filedirs = FigureQuadraticApprox(α,β,x₀,xSol,ε,itmax=itmax,methods=[mtd],mrk_color = [:red], max_iter_plotted=31)
# annotate!([(0.35, .75 ,text(L"\phi(t) = |t|^\alpha - \beta",12))])
# label_points!(pltCARM2,filedirs[1],5,print_proj=true,num_print_proj=2,ProjName="P^S")
# ##

# mtd = :MAP
# pltMAAP2,filedirs = FigureQuadraticApprox(α,β,x₀,xSol,ε,itmax=itmax,methods=[mtd], mrk_color = [:green], max_iter_plotted=15)
# annotate!([(0.35, .75 ,text(L"\phi(t) = |t|^\alpha - \beta",12))])
# label_points!(pltMAAP2,filedirs[1],7,print_proj=true,var_name="z",num_print_proj=2,ProjName="P^S")




# ##
# plt = plot(size=(800,800),pltMAAP1,pltCARM1,pltMAAP2,pltCARM2)
# pltLackEB = plot(size=(800,400),pltMAAP1,pltCARM1)
# pltEB = plot(size=(800,400),pltMAAP2,pltCARM2)
# ## Save Figs
# savefig(pltMAAP1,plotsdir("AABBIS21FigApprox_Quadratic_LackEB_MAP.pdf"))
# savefig(pltCARM1,plotsdir("AABBIS21FigApprox_Quadratic_LackEB_CRM.pdf"))
# savefig(pltMAAP2,plotsdir("AABBIS21FigApprox_Quadratic_EB_MAP.pdf"))
# savefig(pltCARM2,plotsdir("AABBIS21FigApprox_Quadratic_EB_CRM.pdf"))
# savefig(plt,plotsdir("AABBIS21FigApprox_Quadratic_LackEB_EB.pdf"))
# savefig(pltLackEB,plotsdir("AABBIS21FigApprox_Quadratic_LackEB.pdf"))
# savefig(pltEB,plotsdir("AABBIS21FigApprox_Quadratic_EB.pdf"))
# ##
