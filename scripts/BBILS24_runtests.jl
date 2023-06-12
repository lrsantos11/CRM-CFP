__precompile__()
using DrWatson
@quickactivate "CRM-CFP"
include(srcdir("CRM-CFP.jl"))


"""juliadocs

[^Behling2024] R. Behling, J.-Y. Bello-Cruz, A. N. Iusem, D. Liu e L.-R. Santos, scripts/BBILS23_runtests.jl.
"""

##
include(scriptsdir("BBILS24_Ellipsoids.jl"))
