"""
Creates plot with the path of method. 
"""
function method_path!(plt::Plots.Plot,
                    mat::Array; 
                    λ::Float64 = .97,  
                    kwargs...)
    x = @view mat[:,1]
    y = @view mat[:,2]
    num_arrows = length(x)-1
    for index = 1:num_arrows
        x2 = x[index] + λ*(x[index+1] - x[index])
        y2 = y[index] + λ*(y[index+1] - y[index])
        plot!(plt,[x[index],x2],[y[index],y2], lt = :path,  label="", arrow = :closed, ls = :dash; kwargs...)
    end
    return plt
end

function label_points!(plt::Plots.Plot,
                       mat::Vector;
                       num_points::Int = 1,
                       var_name::String = "z^",
                       var_indices::Vector{Int} = Int[],
                       label_size::Int = 8,
                       xshift::Real = 0.0,
                       yshift::Real = 0.0,
                       kwargs...)
    list_pts = @view mat[1:num_points]
    if isempty(var_indices)
         var_indices = collect(0:num_points -1)
    end
    for pto in eachrow(list_pts)
        ptox = pto[1][1]
        ptoy = pto[1][2]
        annotate!(plt,[(ptox+xshift, ptoy+yshift,text(L"%$(var_name){%$(var_indices[pto.indices[1]])}",label_size))];kwargs...)
    end
    return plt
end


function label_points!(plt::Plots.Plot,filedir::String,num_points::Int=5;var_name::String="x", 
                      print_proj::Bool=false,num_print_proj::Int=1,ProjName::String = "P_K")
    xmtd = (readdlm(filedir))[1:num_points,3:end]
    for pto in  eachrow(xmtd)
        if isodd(pto.indices[1])
            annotate!(plt,[(pto[1], -.075 ,text(L"%$(var_name)^{%$(div(pto.indices[1]-1,2))}",12))])
        elseif print_proj && pto.indices[1] <= 2*num_print_proj
            annotate!(plt,[(pto[1]+.05, pto[2]+0.075 ,text(L"%$(ProjName) (%$(var_name)^{%$(div(pto.indices[1]-1,2))})",8))])
        end    
    end
    return plt
end

