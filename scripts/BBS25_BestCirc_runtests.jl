__precompile__()
using DrWatson
@quickactivate "CRM-CFP"
cpu_model = Sys.cpu_info()[1].model
@info cpu_model
if occursin("Intel", cpu_model) || occursin("AMD", cpu_model)
    global islinux = true
    using MKLSparse
    @info "Using MKL and MKLSparse"
else
    global islinux = false
    using AppleAccelerate
    using ThreadedSparseArrays
end
include(srcdir("CRM-CFP.jl"))


"""juliadocs

[^Behling2024] R. Behling, J.-Y. Bello-Cruz,  and L.-R. Santos. Best approximation point using cicumcenter schemes. 2024. In preparation.
"""

##
include(scriptsdir("BBS25_NearestCorrelationMatrix.jl"))
