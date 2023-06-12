__precompile__
using DrWatson
@quickactivate "CRM-CFP"
include(srcdir("CRM-CFP.jl"))
using CSV, DataFrames
using BenchmarkProfiles, LaTeXStrings, Plots
pgfplotsx()
"""juliadocs
This script builds the plots that are presented in  [^Araujo20]

[^Araujo20] G. Araújo, R. Arefidamghani, R. Behling, J.-Y. Bello-Cruz, A. N. Iusem, e L.-R. Santos, 
"""
##

include("AABBIS21_Figs.jl")
include("AABBSI21_ProjEllipsoids.jl")