include("auxplots.jl")
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
pltCRM1 = FigureQuadratic(α,β,x₀,xSol,ε,itmax=itmax,methods=[mtd],mrk_color = [:red], max_iter_plotted=5)
annotate!([(0.35, .75 ,text(L"\phi(t) = |t|^\alpha",12))])
print_points!(pltCRM1,mtd,7,print_proj=true,num_print_proj=2)
##
mtd = :MAP
pltMAP1 = FigureQuadratic(α,β,x₀,xSol,ε,itmax=itmax,methods=[mtd], mrk_color = [:green], max_iter_plotted=21)
annotate!([(0.35, .75 ,text(L"\phi(t) = |t|^\alpha",12))])
print_points!(pltMAP1,mtd,7,print_proj=true,var_name="z",num_print_proj=2)

##
# Quadratic with error bound \phi(t) = αx² + β 

α, β = 1., -0.06
xSol = [sqrt(abs(β)/α), 0.0]
##
mtd = :CRM
pltCRM2 = FigureQuadratic(α,β,x₀,xSol,ε,itmax=itmax,methods=[mtd],mrk_color = [:red], max_iter_plotted=31)
annotate!([(0.35, .75 ,text(L"\phi(t) = |t|^\alpha - \beta",12))])
print_points!(pltCRM2,mtd,5,print_proj=true,num_print_proj=2)
##

mtd = :MAP
pltMAP2 = FigureQuadratic(α,β,x₀,xSol,ε,itmax=itmax,methods=[mtd], mrk_color = [:green], max_iter_plotted=15)
annotate!([(0.35, .75 ,text(L"\phi(t) = |t|^\alpha - \beta",12))])
print_points!(pltMAP2,mtd,7,print_proj=true,var_name="z",num_print_proj=2)




##
plt = plot(size=(800,800),pltMAP1,pltCRM1,pltMAP2,pltCRM2)
pltLackEB = plot(size=(800,400),pltMAP1,pltCRM1)
pltEB = plot(size=(800,400),pltMAP2,pltCRM2)
## Save Figs
savefig(pltMAP1,plotsdir("FigureQuadraticLackEB_MAP.pdf"))
savefig(pltCRM1,plotsdir("FigureQuadraticLackEB_CRM.pdf"))
savefig(pltMAP2,plotsdir("FigureQuadraticEB_MAP.pdf"))
savefig(pltCRM2,plotsdir("FigureQuadraticEB_CRM.pdf"))
savefig(plt,plotsdir("FigureQuadraditcLackEB_EB.pdf"))
savefig(pltLackEB,plotsdir("FigureQuadraticLackEB.pdf"))
savefig(pltEB,plotsdir("FigureQuadraticEB.pdf"))
##
savefig(pltLackEB,"../../../Draft/New/CRM-Convergence-CFP/FigureQuadraticLackEB.pdf")
savefig(pltEB,"../../../Draft/New/CRM-Convergence-CFP/FigureQuadraticEB.pdf")