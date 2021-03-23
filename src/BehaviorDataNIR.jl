module BehaviorDataNIR

using Impute, LinearAlgebra, ProgressMeter, PyPlot, HDF5, Images, Statistics, FlavellBase, VideoIO, HMMBase

include("unit.jl")
include("util.jl")
include("stage_data.jl")
include("cam_data.jl")
include("sync.jl")

export 
    # stage_data.jl
    zero_stage, 
    impute_list,
    speed,
    time_axis,
    Δpos_angle,
    angular_velocity,
    ang_btw_vec,
    reversal_state,
    # unit.jl
    unit_stage_unit_to_mm,
    unit_bfs_pix_to_mm,
    # cam_data.jl
    nmp_vec,
    # util.jl
    diff_lag,
    vec_ang,
    read_h5,
    # sync.jl
    sync_timing,
    sync_stim,
    signal_stack_repeatability

end # module
