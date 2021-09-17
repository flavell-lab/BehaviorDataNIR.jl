module BehaviorDataNIR

using Impute, LinearAlgebra, ProgressMeter, HDF5, Images, Statistics,
    FlavellBase, UNet2D, SegmentationStats, StatsBase,
    Combinatorics, Interpolations, PyCall, Optim, VideoIO, Luxor

include("init.jl")
include("analysis.jl")
include("unit.jl")
include("util.jl")
include("stage_data.jl")
include("cam_data.jl")
include("sync.jl")
include("segmentation.jl")
include("spline_data.jl")
include("video.jl")

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
    get_reversal_events,
    compute_reversal_times,
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
    read_pos_feature,
    read_stage,
    recenter_angle,
    local_recenter_angle,
    vec_to_angle,
    make_vec,
    get_lsqerr,
    savitzky_golay_filter,
    euclidean_dist,
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
    compute_worm_spline!,
    compute_worm_thickness,
    get_segment_end_matrix,
    # spline_data.jl
    get_worm_body_angle,
    get_worm_vector,
    get_tot_worm_curvature,
    self_intersect_ratio,
    # analysis.jl
    get_tuning,
    reversal_neuron_model,
    RIM_model,
    ewma,
    forward_neuron_model,
    turning_neuron_model,
    # video.jl,
    encode_movie,
    write_behavior_video,
    add_text_to_image

end # module
