using Plots
using LaTeXStrings
pgfplotsx()
plotlyjs
"""
Creates plot with the path of method. 
"""
function MethodPath(plt::Plots.Plot,mat::Array;ltype = :dash,lwidth=.5,color=:dodgerblue2)
    x = @view mat[:,1]
    y = @view mat[:,2]
    num_arrows = length(x)-1
    λ = .97
    for index = 1:num_arrows
        x2 = x[index] + λ*(x[index+1] - x[index])
        y2 = y[index] + λ*(y[index+1] - y[index])
        plot!(plt,[x[index],x2],[y[index],y2],st=:path,ls=ltype,c=color,label=false,lw=lwidth,arrow=:closed)
    end
end

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
    for (index,mtd) in enumerate(methods)
        filedir = datadir("sims","x"*String(mtd)*".dat")
        func = getfield(Main,mtd)
        Result = func(x₀, ProjectARight,ProjectBRight,itmax=itmax,filedir=filedir,EPSVAL=ε,xSol=xSol)
        @show Result
        pts_read_mtd = 2*min(Result.iter_total,max_iter_plotted) + 1
        xmtd = (readdlm(filedir))[1:pts_read_mtd,:]
        scatter_points_Proj_x = @view xmtd[2:2:end,1]
        scatter_points_Proj_y = @view xmtd[2:2:end,2]
        scatter_points_mtd_x = @view xmtd[3:2:end,1]
        scatter_points_mtd_y = @view xmtd[3:2:end,2]
        scatter!(scatter_points_mtd_x,scatter_points_mtd_y, c=mrk_color[index], marker=(3,:circle),
                label="$(String(mtd)) ($(Result.iter_total) it.)")
        scatter!(scatter_points_Proj_x,scatter_points_Proj_y, c=:blue, marker=(2,:circle),
                label="")                
        MethodPath(plt,xmtd)
    end
    return plt
end

function print_points!(plt::Plots.Plot,mtd::Symbol,num_points::Int=5;var_name::String="x", 
                      print_proj::Bool=false,num_print_proj::Int=1)
    filedir = datadir("sims","x"*String(mtd)*".dat")
    xmtd = (readdlm(filedir))[1:num_points,:]
    for pto in  eachrow(xmtd)
        if isodd(pto.indices[1])
            annotate!(plt,[(pto[1], -.075 ,text(L"%$(var_name)^{%$(div(pto.indices[1]-1,2))}",12))])
        elseif print_proj && pto.indices[1] <= 2*num_print_proj
            annotate!(plt,[(pto[1]-.05, pto[2]+0.075 ,text(L"P_K(%$(var_name)^{%$(div(pto.indices[1]-1,2))})",8))])
        end    
    end
    return plt
end