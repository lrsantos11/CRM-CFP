using Plots
using LaTeXStrings
pgfplotsx(size=(400,400),fg_colour_axis=:lightgray,framestyle=:origin,aspect_ratio=:equal,draw_arrow=true)

"""
Creates plot with the path of method. 
"""
function MethodPath(plt::Plots.Plot,mat::Array;color=:dodgerblue2,ltype = :dash,lwidth=.5)
    x = mat[1][1:end-1]
    y = mat[2][1:end-1]
    vx = .97(mat[1][2:end] - x)
    vy = .97(mat[2][2:end] - y)
    quiver!(plt,x,y,quiver=(vx,vy),arrow=true,ls=ltype,c=color,lw=lwidth,label="")
end
##

function FigureQuadratic(α::Number,β::Number, x₀::Vector, xSol::Vector,ε::Number, itmax::Int64=6)
##
    # On the Right
    # f(t) = α t^2 + β 
    shiftvec = [0,β]
    yH = 0.0
    Hyper = IndAffine([0,1.], 0.0)
    ProjectARight(x) =  ProjectEpiQuadratic(x-shiftvec,α) + shiftvec
    ProjectBRight(x) =  ProjectIndicator(Hyper,x)
    ReflectARight(x) = Reflection(x,ProjectARight)
    ReflectBRight(x) =  ReflectIndicator(Hyper,x)
    #MAP
    MAPfile = datadir("sims","xMAP.dat")
    MAPResult = MAP(x₀, ProjectARight,ProjectBRight,itmax=itmax,filedir=MAPfile,EPSVAL=ε,xSol=xSol)
    xMAPDir = (readdlm(MAPfile))[1:20,:]
    
    #CRM
    CRMfile = datadir("sims","CRM.dat")
    CRMResult = CRM(x₀, ReflectARight,ReflectBRight,itmax=itmax,filedir=CRMfile,EPSVAL=ε,xSol=xSol)
    xCRMDir = readdlm(CRMfile)
##
    
    # Plots
    ymax = 1.25
    plt = plot(x -> α*x^2 + β,-0.5,sqrt((ymax - β)/α),lw=1,c=:blue,fill = (ymax, 0.5, :dodgerblue),label="")
    scatter!([x₀[1]],[x₀[2]],label="",marker=:diamond)
    MAPscatter = [xMAPDir[:,1],xMAPDir[:,2]]
    scatter!(MAPscatter[1][2:end],MAPscatter[2][2:end], label="MAP ($(MAPResult.iter_total) it.)")
    MethodPath(plt,MAPscatter,color=:green)
    CRMscatter = [xCRMDir[:,1],xCRMDir[:,2]]
    scatter!(CRMscatter[1][2:end],CRMscatter[2][2:end], label="CRM ($(CRMResult.iter_total) it.)",marker=:square)
    MethodPath(plt,CRMscatter,color=:yellow)
    annotate!([(x₀[1], yH -.1 ,text(L"x_0",12))])
    if β  ≈ 0.0
        annotate!([(0.35, .75 ,text(L"f(t) = | t|^\alpha",12))])
    else
        annotate!([(0.35, .75 ,text(L"f(t) = | t|^\alpha + \beta",12))])
    end
    ##
    # annotate!([(-0.4, .9 ,text(L"X",12))])
    # plot!(x->yH,-2.2,1.2,c=:black,lw=1,fill = (-.64+yH, 0.8, :ivory3))
    # plot!(x->-(x+2.2)^2 + yH,-3.,-2.2,c=:black,lw=1,fill = (-.64+yH, 0.8, :ivory3))
    # plot!(x->-2(x-1.2)^2 + yH,1.2,sqrt(.32)+1.2,c=:black,lw=1,fill = (-.64+yH, 0.8, :ivory3))
    # plot!(x->yH,-0.5,0.,c=:red,lw=2);
    # plot!(x->0,-0.5,0.,c=:red,lw=2,fill = (1.5, 0.5, :dodgerblue));
    # scatter!([1,-2],[.5+yH,.5+yH],c=:white,lc=:black,lw=2,marker=(3,:circle))
   
    return plt
end
