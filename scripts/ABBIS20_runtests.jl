__precompile__()
using DrWatson
@quickactivate "CRM-CFP"
include(srcdir("CRM-CFP.jl"))
using CSV, DataFrames
using BenchmarkProfiles, LaTeXStrings, Plots
pgfplotsx()

"""juliadocs
This script builds the plots that are presented in Figure 4.1 of Example 4.12 from [^Arefidamghani20]

[^Arefidamghani20] R. Arefidamghani, R. Behling, J.-Y. Bello-Cruz, A. N. Iusem, e L.-R. Santos, 
“The circumcentered-reflection method achieves better rates than alternating projections”, 
arXiv, ago. 2020, [Online]. Disponível em: http://arxiv.org/abs/2007.14466.
"""


"""
    FigureQuadratic(α,β,x₀,xSol,ε)
"""
function FigureQuadratic(α::Number,β::Number, x₀::Vector, xSol::Vector,ε::Number; itmax::Int64=6,  
    methods = [:MAP, :CRM], mrk_color = [:green, :blue], max_iter_plotted ::Int=13, ymax::Number = 1.25, framestyle=:origin)
    # On the Right
    # f(t) = α t^2 + β 
    shiftvec = [0,β]
    Hyper = IndAffine([0,1.], 0.0)
    ProjectARight(x) =  ProjectEpiQuadratic(x-shiftvec,α) + shiftvec
    ProjectBRight(x) =  ProjectIndicator(Hyper,x)
    plt = plot(size=(400,400), fg_colour_axis=:lightgray,
              framestyle=framestyle,
               aspect_ratio=:equal,draw_arrow=true,ticks=:none,grid=:none,legend=:topright)
    plot!(x -> α*x^2 + β,-0.5,sqrt((ymax - β)/α),lw=1,c=:blue,fill = (ymax, 0.5, :dodgerblue),label="")
    scatter!([x₀[1]],[x₀[2]],label="",marker=(3,:circle))
    timenow = Dates.now()
    filedirs = String[]
    for (index,mtd) in enumerate(methods)
        filename = savename("ABBIS20Fig4-1",(mtd=mtd,date=timenow),"csv",sort=false)
        filedir = datadir("sims",filename)
        push!(filedirs,filedir)
        func = getfield(Main,mtd)
        Result = func(x₀, ProjectARight,ProjectBRight,itmax=itmax,filedir=filedir,EPSVAL=ε,xSol=xSol,print_intermediate = true)
        @show Result
        pts_read_mtd = min(Result.iter_total,max_iter_plotted) + 1
        xmtd = (readdlm(filedir))[1:pts_read_mtd,3:end]
        scatter_points_Proj_x = @view xmtd[2:2:end,1]
        scatter_points_Proj_y = @view xmtd[2:2:end,2]
        scatter_points_mtd_x = @view xmtd[3:2:end,1]
        scatter_points_mtd_y = @view xmtd[3:2:end,2]
        scatter!(scatter_points_Proj_x,scatter_points_Proj_y, c=:blue, marker=(2,:circle),
                label="")                
        scatter!(scatter_points_mtd_x,scatter_points_mtd_y, c=mrk_color[index], marker=(3,:circle),
                label="$(String(mtd)) ($(Result.iter_total) it.)")        
        MethodPath!(plt, xmtd, color=:dodgerblue, alpha=0.2)
    end
    return plt,filedirs
end

##
#
itmax = 10^6
ε = 1e-3
x₀ = [1.1,0]
##
# Quadratic without error bound
α, β = 1., 0.0
xSol = [0.0, 0.0]
# 
mtd = :CRM
pltCRM1,filedirs = FigureQuadratic(α,β,x₀,xSol,ε,itmax=itmax,methods=[mtd],mrk_color = [:red], max_iter_plotted=5)

annotate!([(0.35, .75 ,text(L"\phi(t) = |t|^\alpha",12))])
label_points!(pltCRM1, filedirs[1], 7, print_proj=true, num_print_proj=2, xProjshift=-0.05)

##
mtd = :MAP
pltMAP1,filedirs = FigureQuadratic(α,β,x₀,xSol,ε,itmax=itmax,methods=[mtd], mrk_color = [:green], max_iter_plotted=21)
annotate!([(0.35, .75 ,text(L"\phi(t) = |t|^\alpha",12))])
label_points!(pltMAP1,filedirs[1],7,print_proj=true,var_name="z",num_print_proj=2)

##
# Quadratic with error bound \phi(t) = αx² + β 

α, β = 1., -0.06
ymax = 1.25
max_iter_plotted = 5
xSol = [sqrt(abs(β)/α), 0.0]
##
mtd = :CRM
pltCRM2, filedirs = FigureQuadratic(α, β, x₀, xSol, ε, itmax=itmax, methods=[mtd], mrk_color=[:red], max_iter_plotted = max_iter_plotted, framestyle = :none)
hline!([0.0], lw=2, c=:black, label="")
annotate!([(0., .55 ,text(L"K",12))])
annotate!([(-0.4, -0.075, text(L"U", 12))])
label_points!(pltCRM2,filedirs[1],5, print_proj=true,num_print_proj = 2, xProjshift = -0.06)
savefig(pltCRM2, plotsdir("Iusem_QuadraticEB_CRM_1.pdf"))

##
Proj_x₀ = (readdlm(filedir))[2,3:4]
ReflecK_x₀ = 2*Proj_x₀ - x₀
Hyper = IndAffine([0, 1.0], 0.0)
ReflecU_x₀ = ReflectIndicator(Hyper, ReflecK_x₀)

scatter!([Singleton(v) for v in [ReflecK_x₀, ReflecU_x₀]], c=:blue, marker=(2, :circle), label="")
annotate!([(ReflecK_x₀[1] - 0.06, ReflecK_x₀[2] + 0.075, text(L"R_K(x^0)", 8))])
annotate!(ReflecU_x₀[1] - 0.2, ReflecU_x₀[2] + 0.05, text(L"R_UR_K(x^0)", 8))
MethodPath!(pltCRM2, [x₀'; ReflecK_x₀'; ReflecU_x₀';  x₀'], color=:dodgerblue, alpha=0.2)
savefig(pltCRM2, plotsdir("Iusem_QuadraticEB_CRM_2.pdf"))

##

mtd = :MAP
pltMAP2,filedirs = FigureQuadratic(α,β,x₀,xSol,ε,itmax=itmax,methods=[mtd], mrk_color = [:green], max_iter_plotted=15, framestyle = :none)
annotate!([(0.0, 0.55, text(L"K", 12))])
annotate!([(-0.4, -0.075, text(L"U", 12))])
hline!([0.0], lw=2, c=:black, label="")
label_points!(pltMAP2, filedirs[1], 7, print_proj=true, var_name="z", num_print_proj=2, xProjshift=-0.06)
savefig(pltMAP2, plotsdir("Iusem_QuadraticEB_MAP_1.pdf"))

##
scatter!([Singleton(v) for v in [ReflecK_x₀, ReflecU_x₀]], c=:blue, marker=(2, :circle), label="")
annotate!([(ReflecK_x₀[1] - 0.06, ReflecK_x₀[2] + 0.075, text(L"R_K(z^0)", 8))])
annotate!(ReflecU_x₀[1] - 0.2, ReflecU_x₀[2] + 0.05, text(L"R_UR_K(z^0)", 8))
MethodPath!(pltMAP2, [x₀'; ReflecK_x₀'; ReflecU_x₀'; x₀'], color=:dodgerblue, alpha=0.2)
savefig(pltMAP2, plotsdir("Iusem_QuadraticEB_MAP_2.pdf"))



##
plt = plot(size=(800,800),pltMAP1,pltCRM1,pltMAP2,pltCRM2)
pltLackEB = plot(size=(800,400),pltMAP1,pltCRM1)
pltEB = plot(size=(800,400),pltMAP2,pltCRM2)
## Save Figs
savefig(pltMAP1,plotsdir("ABBIS20Fig4-1_QuadraticLackEB_MAP.pdf"))
savefig(pltCRM1,plotsdir("ABBIS20Fig4-1_QuadraticLackEB_CRM.pdf"))
savefig(pltMAP2,plotsdir("ABBIS20Fig4-1_QuadraticEB_MAP.pdf"))
savefig(pltCRM2,plotsdir("ABBIS20Fig4-1_QuadraticEB_CRM.pdf"))
savefig(plt,plotsdir("ABBIS20Fig4-1_QuadraditcLackEB_EB.pdf"))
savefig(pltLackEB,plotsdir("ABBIS20Fig4-1_QuadraticLackEB.pdf"))
savefig(pltEB,plotsdir("ABBIS20Fig4-1_QuadraticEB.pdf"))
##
