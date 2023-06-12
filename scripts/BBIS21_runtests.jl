__precompile__
using DrWatson
@quickactivate "CRM-CFP"
include(srcdir("CRM-CFP.jl"))
using CSV, DataFrames
using BenchmarkProfiles,  LazySets
pgfplotsx()
# include(scriptsdir("plots_util.jl"))

"""juliadocs

[^Behling2023] R. Behling, J.-Y. Bello-Cruz, A. N. Iusem, e L.-R. Santos, On the centralization of the circumcentered-reflection method.
"""
##

# include(scriptsdir("BBIS21_Figs.jl"))
##
# include(scriptsdir("BBIS21_Ellipsoids.jl"))


# include(scriptsdir("BBIS21_Sec4.2.jl"))