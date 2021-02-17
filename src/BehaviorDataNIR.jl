module BehaviorDataNIR

using Impute, LinearAlgebra, PyCall, PyPlot, HDF5, Images, Statistics, ProgressMeter,  FlavellBase, VideoIO

ski_morphology = pyimport_conda("skimage.morphology", "scikit-image")
np = pyimport("numpy")
cv2 = pyimport("cv2");

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
    standardize,
    reshape_array,
    eval_unet,
    create_model,
    unet,
    overlay_images,
    save_images,
    encode_movie

end # module
