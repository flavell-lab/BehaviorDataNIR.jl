module BehaviorDataNIR

using Impute, LinearAlgebra, ProgressMeter, PyPlot, HDF5, Images, Statistics,  FlavellBase, VideoIO, HMMBase

include("unit.jl")
include("util.jl")
include("stage_data.jl")
include("cam_data.jl")


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
    make_vec,
# unit.jl
    unit_stage_unit_to_mm,
    unit_bfs_pix_to_mm,
# cam_data.jl
    read_hdf5, 
    nmp_vec,
# util.jl
    diff_lag,
    vec_ang

end # module
