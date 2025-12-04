__precompile__()

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
    @info "Using AppleAccelerate and ThreadedSparseArrays"
end


using DrWatson

const global ZERO_VAL = 1e-15
include("CRM_utils.jl")
include("Plots_util.jl")
# include("read_netlib_MPS.jl")
include("Ellipsoids_utils.jl")
include("MAP.jl")
include("DRM.jl")
include("SPM.jl")
include("CRM.jl")
include("CentCRM.jl")
include("SucCentCRM.jl")
include("PACA.jl")
include("paralellCRM.jl")