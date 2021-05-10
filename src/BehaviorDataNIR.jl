module BehaviorDataNIR

using Impute, LinearAlgebra, ProgressMeter, HDF5, Images, Statistics,
    FlavellBase, HMMBase, UNet2D, SegmentationStats,
    Combinatorics, LinearAlgebra, Interpolations, PyCall, Optim

include("unit.jl")
include("util.jl")
include("stage_data.jl")
include("cam_data.jl")
include("sync.jl")
include("segmentation.jl")
include("spline_data.jl")
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
    offset_xy,
    # unit.jl
    unit_bfs_pix_to_stage_unit,
    unit_stage_unit_to_mm,
    unit_bfs_pix_to_mm,
    # cam_data.jl
    nmp_vec,
    # util.jl
    diff_lag,
    vec_ang,
    read_h5,
    recenter_angle,
    vec_to_angle,
    make_vec,
    savitzky_golay_filter,
    # sync.jl
    sync_timing,
    sync_stim,
    signal_stack_repeatability,
    nir_vec_to_confocal,
    unlag_vec,
    nir_to_confocal_t,
    # segmentation.jl
    downsample_unet_input,
    segment_worm!,
    medial_axis,
    fit_spline,
    # spline_data.jl
    get_worm_body_angle,
    get_worm_vector,
    # analysis.jl
    get_tuning

end # module
