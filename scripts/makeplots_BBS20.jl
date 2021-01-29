"""
This script builds the plot that is presented in Figure 1 in [^Behling2020]

[^Behling2020] R. Behling, J.-Y. Bello-Cruz, e L.-R. Santos, “On the Circumcentered-Reflection Method 
for the Convex Feasibility Problem”, Numer. Algorithms, jul. 2020, doi: [10.1007/s11075-020-00941-6](https://doi.org/10.1007/s11075-020-00941-6). 
"""

include("plots_util.jl")




function plotBalls(x₀::Vector,ProjectA::Function,ProjectB::Function; 
            xSol::Vector=[], ε::Number=1e-6, itmax::Int64=10, methods = [:MAP, :CRM],
            mrk_color = [:green, :blue], mrk_type = [:square, :diamond])
    timenow = Dates.now()
    filedirs = String[]
    for (index,mtd) in enumerate(methods)
        filename = savename("BBS20Fig1",(mtd=mtd,date=timenow),"csv",sort=false)
        filedir = datadir("sims",filename)
        push!(filedirs,filedir)
        func = getfield(Main,mtd)
        Result = func(x₀, ProjectA,ProjectB,itmax=itmax,filedir=filedir,EPSVAL=ε,xSol=xSol)
        @show Result
        xmtd = readdlm(filedir)
        scatter_points_mtd_x = @view xmtd[:,1]
        scatter_points_mtd_y = @view xmtd[:,2]
        scatter!(scatter_points_mtd_x,scatter_points_mtd_y, c=mrk_color[index], marker=(4,mrk_type[index]), 
                line=(:path,.7), label="$(String(mtd)) ($(Result.iter_total) it.)") 
        # scatter!(scatter_points_mtd_x,scatter_points_mtd_y, c=mrk_color[index], marker=(3,:circle),
                # label="$(String(mtd)) ($(Result.iter_total) it.)")   
    end
    return plt,filedirs
end


v1 = [.8,0]
r1 = 1.
v2 = [-.8,0]
r2 = 1.
plt = plot(size=(400,400),fg_colour_axis=:lightgray,framestyle=:none,
               aspect_ratio=:equal,draw_arrow=true,ticks=:none,grid=:none)
 plot!(circleShape(v1,r1), seriestype=:shape, lw = 1, c=:blue, linecolor=:blue,  
    fillalpha = .4, label="")
plot!(circleShape(v2,r2), seriestype=:shape, lw = 1, c=:green, linecolor=:green,  
    fillalpha = .4, label="")
annotate!([(v1[1],  0,text(L"X_2",16,:bottom, color=:blue))])
annotate!([(v2[1],  0,text(L"X_1",16,:bottom, color=:green))])
    
ProjectA(x) =  ProjectBall(x,v2, r2)
ProjectB(x) =  ProjectBall(x,v1, r1)

τ = .3
xzero = Float64[2.1,-2.8]
plotBalls(xzero,ProjectA,ProjectB)
annotate!([(xzero[1],  xzero[2]-τ,text(L"x^0",14))])
ylims!(xzero[2]-2*τ,max(v1[2],v2[2]) + max(r1,r2) +τ)

savefig(plt,plotsdir("BBS20Fig1_L2Balls.pdf"))
