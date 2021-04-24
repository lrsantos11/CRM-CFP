using DrWatson
@quickactivate "CRM-CFP"
include(srcdir("CRM-CFP.jl"))

using DataFrames
using LaTeXStrings, Plots
pgfplotsx()


# include("runtests_ABBIS20.jl")
# include("runtests_BBS20.jl")
include("runtests_AABBIS20.jl")



