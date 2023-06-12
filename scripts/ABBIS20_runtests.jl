__precompile__()
using DrWatson
@quickactivate "CRM-CFP"
include(srcdir("CRM-CFP.jl"))
using CSV, DataFrames
using BenchmarkProfiles, LaTeXStrings, Plots
# pgfplotsx()

"""juliadocs
This script builds the plots that are presented in Figure 4.1 of Example 4.12 from [^Arefidamghani20]

[^Arefidamghani20] R. Arefidamghani, R. Behling, J.-Y. Bello-Cruz, A. N. Iusem, e L.-R. Santos, 
“The circumcentered-reflection method achieves better rates than alternating projections”, 
arXiv, ago. 2020, [Online]. Disponível em: http://arxiv.org/abs/2007.14466.
"""

include("plots_util.jl")

"""
    FigureQuadratic(α,β,x₀,xSol,ε)
"""
function FigureQuadratic(α::Number,β::Number, x₀::Vector, xSol::Vector,ε::Number; itmax::Int64=6,  
    methods = [:MAP, :CRM], mrk_color = [:green, :blue], max_iter_plotted ::Int=13, ymax::Number = 1.25)
    # On the Right
    # f(t) = α t^2 + β 
    shiftvec = [0,β]
    Hyper = IndAffine([0,1.], 0.0)
    ProjectARight(x) =  ProjectEpiQuadratic(x-shiftvec,α) + shiftvec
    ProjectBRight(x) =  ProjectIndicator(Hyper,x)
    plt = plot(size=(400,400),fg_colour_axis=:lightgray,framestyle=:origin,
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
        pts_read_mtd = 2*min(Result.iter_total,max_iter_plotted) + 1
        xmtd = (readdlm(filedir))[1:pts_read_mtd,3:end]
        scatter_points_Proj_x = @view xmtd[2:2:end,1]
        scatter_points_Proj_y = @view xmtd[2:2:end,2]
        scatter_points_mtd_x = @view xmtd[3:2:end,1]
        scatter_points_mtd_y = @view xmtd[3:2:end,2]
        scatter!(scatter_points_Proj_x,scatter_points_Proj_y, c=:blue, marker=(2,:circle),
                label="")                
        scatter!(scatter_points_mtd_x,scatter_points_mtd_y, c=mrk_color[index], marker=(3,:circle),
                label="$(String(mtd)) ($(Result.iter_total) it.)")        
       method_path(plt,xmtd)
    end
    return plt,filedirs
end


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
label_points!(pltCRM1,filedirs[1],7,print_proj=true,num_print_proj=2)
##
mtd = :MAP
pltMAP1,filedirs = FigureQuadratic(α,β,x₀,xSol,ε,itmax=itmax,methods=[mtd], mrk_color = [:green], max_iter_plotted=21)
annotate!([(0.35, .75 ,text(L"\phi(t) = |t|^\alpha",12))])
label_points!(pltMAP1,filedirs[1],7,print_proj=true,var_name="z",num_print_proj=2)

##
# Quadratic with error bound \phi(t) = αx² + β 

α, β = 1., -0.06
xSol = [sqrt(abs(β)/α), 0.0]
##
mtd = :CRM
pltCRM2,filedirs = FigureQuadratic(α,β,x₀,xSol,ε,itmax=itmax,methods=[mtd],mrk_color = [:red], max_iter_plotted=31)
annotate!([(0.35, .75 ,text(L"\phi(t) = |t|^\alpha - \beta",12))])
label_points!(pltCRM2,filedirs[1],5,print_proj=true,num_print_proj=2)
##

mtd = :MAP
pltMAP2,filedirs = FigureQuadratic(α,β,x₀,xSol,ε,itmax=itmax,methods=[mtd], mrk_color = [:green], max_iter_plotted=15)
annotate!([(0.35, .75 ,text(L"\phi(t) = |t|^\alpha - \beta",12))])
label_points!(pltMAP2,filedirs[1],7,print_proj=true,var_name="z",num_print_proj=2)




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
