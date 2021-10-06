__precompile__()
include("crm_utils.jl")
global const ZERO_VAL = 1e-15
include("MAP.jl")
include("DRM.jl")
include("SPM.jl")
include("CRM.jl")
include("CentCRM.jl")