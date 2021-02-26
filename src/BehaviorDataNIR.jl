module BehaviorDataNIR

using Impute, LinearAlgebra, PyCall, PyPlot, HDF5, Images, Statistics, ProgressMeter,  FlavellBase, VideoIO

include("unit.jl")
include("util.jl")
include("stage_data.jl")
include("cam_data.jl")
include("worm_segmentation.jl") 

# stage_data.jl
export zero_stage,
    impute_list,
    impute_stage,
    speed,
    time_axis,
    Δpos_angle,
    angular_velocity,
    ang_btw_vec,
    reverse_vec,
    mov_vec,
    cluster,
# unit.jl
    unit_stage_unit_to_mm,
    unit_bfs_pix_to_mm,
# cam_data.jl
    read_hdf5,
    nmp_vec

# TODO: move to UNet2D library
# ski_morphology = pyimport_conda("skimage.morphology", "scikit-image")
# np = pyimport("numpy")
# cv2 = pyimport("cv2");

# # worm_segmentation.jl
#     standardize,
#     reshape_array,
#     eval_unet,
#     create_model,
#     unet,
#     overlay_images,
#     save_images,
#     encode_movie

end # module
