module BehaviorDataNIR

using Impute, LinearAlgebra, ProgressMeter, HDF5, Images, Statistics,
    FlavellBase, HMMBase, UNet2D, SegmentationStats,
    Combinatorics, LinearAlgebra, Dierckx, PyCall

include("unit.jl")
include("util.jl")
include("stage_data.jl")
include("cam_data.jl")
include("sync.jl")
include("segmentation.jl")
include("init.jl")

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
    signal_stack_repeatability,
    # segmentation.jl
    downsample_unet_input,
    segment_worm!,
    medial_axis,
    fit_spline

end # module
