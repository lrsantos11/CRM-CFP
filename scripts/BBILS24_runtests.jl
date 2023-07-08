__precompile__()
using DrWatson
@quickactivate "CRM-CFP"
if occursin("Intel", Sys.cpu_info()[1].model)
    using MKL
end
include(srcdir("CRM-CFP.jl"))


"""juliadocs

[^Behling2024] R. Behling, J.-Y. Bello-Cruz, A. N. Iusem, D. Liu e L.-R. Santos.
"""

##
include(scriptsdir("BBILS24_Ellipsoids.jl"))
