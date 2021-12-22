__precompile__
using DrWatson
@quickactivate "CRM-CFP"
include(srcdir("CRM-CFP.jl"))
using CSV, DataFrames
using BenchmarkProfiles, LaTeXStrings, Plots, LazySets
pgfplotsx()
include(scriptsdir("plots_util.jl"))

"""juliadocs

[^BBIS21] R. Behling, J.-Y. Bello-Cruz, A. N. Iusem, e L.-R. Santos, On the centralization of the circumcentered-reflection method.
"""
##

include(scriptsdir("BBIS21_Figs.jl"))
include(scriptsdir("BBIS21_Ellipsoids.jl"))