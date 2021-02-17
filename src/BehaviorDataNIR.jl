module BehaviorDataNIR

using Impute, LinearAlgebra

include("unit.jl")
include("util.jl")
include("stage_data.jl")
include("cam_data.jl")
include("worm_segmentation.jl") 

# stage_data.jl
export zero_stage,
    impute_stage,
    speed,
    Δpos_angle,
    angular_velocity,
# unit.jl
    unit_stage_unit_to_mm,
    unit_bfs_pix_to_mm,
# cam_data.jl
    nmp_vec,
# worm_segmentation.jl
    standardize

end # module
