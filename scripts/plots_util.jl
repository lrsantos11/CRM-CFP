using Plots
using LaTeXStrings
pgfplotsx()

"""
Creates plot with the path of method. 
"""
function MethodPath(plt::Plots.Plot,mat::Array;
    ltype = :dash,lwidth=.5,color=:dodgerblue2,alpha::Float64=1.0, λ::Float64 = .97)
    x = @view mat[:,1]
    y = @view mat[:,2]
    num_arrows = length(x)-1
    for index = 1:num_arrows
        x2 = x[index] + λ*(x[index+1] - x[index])
        y2 = y[index] + λ*(y[index+1] - y[index])
        plot!(plt,[x[index],x2],[y[index],y2],line=(:path,alpha),ls=ltype,c=color,label="",lw=lwidth,arrow=:closed)
    end
end

function label_points!(plt::Plots.Plot,filedir::String,num_points::Int=5;var_name::String="x", 
                      print_proj::Bool=false,num_print_proj::Int=1)
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

"""
circleShape(h::Float64,k::Float64,r::Float64)

Returns the parametric function for a circle with center `v=(h,k)` and radius `r`.
"""
function circleShape(v,r::Float64)
    h = v[1]
    k = v[2]
    θ = LinRange(0.,2*π,500)
    return h .+ r*sin.(θ), k .+ r*cos.(θ)
end