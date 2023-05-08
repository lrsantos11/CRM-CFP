__precompile__()
using DrWatson
@quickactivate "CRM-CFP"
include(srcdir("CRM-CFP.jl"))


"""juliadocs

[^Behling2023a] R. Behling, J.-Y. Bello-Cruz, A. N. Iusem, D. Liu e L.-R. Santos, A successive centralized circumcentered-reflection method for the convex feasibility problem.
"""

##
include(scriptsdir("BBILS23_Ellipsoids.jl"))
