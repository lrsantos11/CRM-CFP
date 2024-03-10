__precompile__()
using DrWatson
@quickactivate "CRM-CFP"
if occursin("Intel", Sys.cpu_info()[1].model)
    using MKL
end
include(srcdir("CRM-CFP.jl"))


"""juliadocs

[^Behling2024] R. Behling, J.-Y. Bello-Cruz,  and L.-R. Santos. Best approximation point using cicumcenter schemes. 2024. In preparation.
"""

##
include(scriptsdir("BBS25_NearestCorrelationMatrix.jl"))
