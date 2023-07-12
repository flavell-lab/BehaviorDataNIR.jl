var documenterSearchIndex = {"docs":
[{"location":"sync/#Confocal-to-NIR-timesyncing-API","page":"Confocal to NIR timesyncing API","title":"Confocal to NIR timesyncing API","text":"","category":"section"},{"location":"sync/","page":"Confocal to NIR timesyncing API","title":"Confocal to NIR timesyncing API","text":"get_timing_info!\nsync_timing\nsync_stim\nnir_vec_to_confocal\nnir_to_confocal_t\nget_timestamps","category":"page"},{"location":"sync/#BehaviorDataNIR.get_timing_info!","page":"Confocal to NIR timesyncing API","title":"BehaviorDataNIR.get_timing_info!","text":"get_timing_info!(data_dict::Dict, param::Dict, path_h5::String, h5_confocal_time_lag::Integer)\n\nInitializes all timing and syncing variables into data_dict::Dict given param::Dict, path_h5::String and the h5_confocal_time_lag::Integer\n\n\n\n\n\n","category":"function"},{"location":"sync/#BehaviorDataNIR.sync_timing","page":"Confocal to NIR timesyncing API","title":"BehaviorDataNIR.sync_timing","text":"sync_timing(path_h5, n_rec=1)\n\nSync NIR and confocal timing. Returns (n_img_nir, daqmx_ai, daqmx_di, img_metadata)\n\nArguments\n\n\n\npath_h5::String: path to h5 file\nn_rec::Int: number of recordings in the file\n\n\n\n\n\n\n\n","category":"function"},{"location":"sync/#BehaviorDataNIR.sync_stim","page":"Confocal to NIR timesyncing API","title":"BehaviorDataNIR.sync_stim","text":"sync_stim(stim, timing_stack, timing_nir)\n\nAlign stim to confocal and NIR timing. Returns (stim_to_confocal, stim_to_nir)\n\nArguments\n\n\n\nstim::Vector: stim signal\ntiming_stack: confocal timing\ntiming_nir: NIR timing\n\n\n\n\n\n\n\n","category":"function"},{"location":"sync/#BehaviorDataNIR.nir_vec_to_confocal","page":"Confocal to NIR timesyncing API","title":"BehaviorDataNIR.nir_vec_to_confocal","text":"nir_vec_to_confocal(vec, confocal_to_nir, confocal_len)\n\nBins NIR behavioral data to match the confocal time points.\n\nArguments:\n\nvec: behavioral data vector. Can be 1D or 2D; if 2D time should be on the columns\nconfocal_to_nir: confocal to NIR time sync\nconfocal_len: length of confocal dataset\n\n\n\n\n\n","category":"function"},{"location":"sync/#BehaviorDataNIR.nir_to_confocal_t","page":"Confocal to NIR timesyncing API","title":"BehaviorDataNIR.nir_to_confocal_t","text":"nir_to_confocal_t(t, nir_to_confocal)\n\nConverts NIR time point t to confocal time point using nir_to_confocal timesync variable.\n\n\n\n\n\n","category":"function"},{"location":"sync/#BehaviorDataNIR.get_timestamps","page":"Confocal to NIR timesyncing API","title":"BehaviorDataNIR.get_timestamps","text":"get_timestamps(path_h5)\n\nGets NIR timestamps from the NIR data file\n\n\n\n\n\n","category":"function"},{"location":"feeding/#Feeding-import-API","page":"Feeding import API","title":"Feeding import API","text":"","category":"section"},{"location":"feeding/","page":"Feeding import API","title":"Feeding import API","text":"import_pumping!","category":"page"},{"location":"feeding/#BehaviorDataNIR.import_pumping!","page":"Feeding import API","title":"BehaviorDataNIR.import_pumping!","text":"import_pumping!(combined_data_dict::Dict, param::Dict, paths_pumping; prefix::String=\"\", is_split=false)\n\nImport pumping data into a combined dataset given param from csv files. Can add a prefix (default empty string) to all confocal variables. If prefix is added, will output a list of pumping for each dataset rather than merging them together.\n\n\n\n\n\n","category":"function"},{"location":"data/#Data-IO-API","page":"Data IO API","title":"Data IO API","text":"","category":"section"},{"location":"data/#Read-data-from-NIR-video-and-tracking-network","page":"Data IO API","title":"Read data from NIR video and tracking network","text":"","category":"section"},{"location":"data/","page":"Data IO API","title":"Data IO API","text":"nmp_vec\nread_h5\nread_pos_feature","category":"page"},{"location":"data/#BehaviorDataNIR.nmp_vec","page":"Data IO API","title":"BehaviorDataNIR.nmp_vec","text":"nmp_vec(pos_feature)\n\nGiven a 3x2xN array pos_feature of 2D coordinates of the nose, metacorpus and pharynx, this function returns three arrays:\n\nmn: an array of shape 2xN representing the vector from the metacorpus to the nose for each frame.\nmp: an array of shape 2xN representing the vector from the metacorpus to the pharynx for each frame.\nmp⊥: an array of shape 2xN representing the vector perpendicular to mp for each frame.\n\n\n\n\n\n","category":"function"},{"location":"data/#Write-video","page":"Data IO API","title":"Write video","text":"","category":"section"},{"location":"data/","page":"Data IO API","title":"Data IO API","text":"encode_movie\nwrite_behavior_video\nadd_text_to_image\nwrite_mip_video","category":"page"},{"location":"data/#BehaviorDataNIR.write_behavior_video","page":"Data IO API","title":"BehaviorDataNIR.write_behavior_video","text":"write_behavior_video(path_h5, path_vid=nothing; fps=20, encoder_options=nothing, downsample=true)\n\nWrites a video of the behavior data Arguments ––––-\n\npath_h5: path of the HDF5 file  \npath_vid: path of the video to be generated. if nothing (default), automatically generates it  \nfps: frame rate of the video (default: 20 to match the acquisition rate)  \ndownsample: if true downsample by factor of 2   \nencoder_options: in named tuple (e.g. (crf=\"23\", preset=\"slow\"))\nvars: variables to display in the video represented as a tuple (varname, value_arr, color). Default nothing which does not display any text\ntext_pos: position of variables text (default (5,5))\ntext_size: size of variables text (default 20)\ntext_font: font of variables text (default Futura)\ntext_spacing: vertical spacing of variables text (default 1.05)\n\nNotes on the encoder options\n\npreset: possible options are: ultrafast, superfast, veryfast, faster, fast, medium – default preset, slow, slower, veryslow.   crf: The range of the CRF scale is 0–51, where 0 is lossless, 23 is the default, and 51 is worst quality possible.   A lower value generally leads to higher quality, and a subjectively sane range is 17–28.   Consider 17 or 18 to be visually lossless or nearly so; it should look the same or nearly the same as the input but it isn't technically lossless.  \n\n\n\n\n\n","category":"function"},{"location":"data/#BehaviorDataNIR.add_text_to_image","page":"Data IO API","title":"BehaviorDataNIR.add_text_to_image","text":"add_text_to_image(img, text_arr, position, colors, fs, font, spacing)\n\nAdds text to an image.\n\nArguments\n\nimg: image to which text is to be added\ntext_arr: array of strings containing the text to be added\nposition: position of the text in the image\ncolors: array of colors for each text string\nfs: font size\nfont: font face\nspacing: vertical spacing between text strings\n\nReturns\n\nnew_img: image with text added\n\n\n\n\n\n","category":"function"},{"location":"data/#BehaviorDataNIR.write_mip_video","page":"Data IO API","title":"BehaviorDataNIR.write_mip_video","text":"write_mip_video(param_path, num_timepts, ch, path_vid=nothing; fps=20, encoder_options=nothing)\n\nWrites maximum-intensity projection video of MHD data.\n\nArguments:\n\nparam_path::Dict: Parameter path dictionary containing path_dir_mhd and get_basename keys\nnum_timepts: Number of timepoints to create video of\nch: Confocal channel to use\npath_vid: Path to output video\nfps (default 20): frames per second\nencoder_options: in named tuple (e.g. (crf=\"23\", preset=\"slow\"))\n\n\n\n\n\n","category":"function"},{"location":"data/#Merge-data-from-multiple-datasets-with-the-same-animal","page":"Data IO API","title":"Merge data from multiple datasets with the same animal","text":"","category":"section"},{"location":"data/","page":"Data IO API","title":"Data IO API","text":"merge_nir_data!","category":"page"},{"location":"data/#BehaviorDataNIR.merge_nir_data!","page":"Data IO API","title":"BehaviorDataNIR.merge_nir_data!","text":"merge_nir_data!(combined_data_dict::Dict, data_dict::Dict, data_dict_2::Dict, param::Dict)\n\nModifies combined_data_dict to add merged NIR data from two datasets data_dict and data_dict_2, together with a parameter file param for the primary dataset.\n\n\n\n\n\n","category":"function"},{"location":"locomotion/#Locomotion-API","page":"Locomotion API","title":"Locomotion API","text":"","category":"section"},{"location":"locomotion/#Velocity-computation","page":"Locomotion API","title":"Velocity computation","text":"","category":"section"},{"location":"locomotion/","page":"Locomotion API","title":"Locomotion API","text":"get_velocity!\nspeed\nreversal_state\nget_reversal_events\ncompute_reversal_times","category":"page"},{"location":"locomotion/#BehaviorDataNIR.get_velocity!","page":"Locomotion API","title":"BehaviorDataNIR.get_velocity!","text":"get_velocity!(data_dict::Dict, param::Dict; prefix::String=\"\")\n\nGets velocity, speed, reversal, and related variables from data_dict and param. Can add a prefix (default empty string) to all confocal variables.\n\n\n\n\n\n","category":"function"},{"location":"locomotion/#BehaviorDataNIR.speed","page":"Locomotion API","title":"BehaviorDataNIR.speed","text":"speed(Δx::Array{<:AbstractFloat,1}, Δy::Array{<:AbstractFloat,1}, Δt::AbstractFloat)\n\nComputes speed in mm/s \n\nArguments\n\nΔx: discrete difference of x\nΔy: discrete difference of y\nΔt: time interval\n\n\n\n\n\nspeed(x, y; lag::Int, fps=FLIR_FPS)\n\nComputes speed using x, y coordinates\n\nlag determines the interval at which to compute the discrete difference\n\nArguments\n\nx: list of x\ny: list of y\nlag: number of time points for discerete difference interval\nfps: fps\n\n\n\n\n\n","category":"function"},{"location":"locomotion/#BehaviorDataNIR.get_reversal_events","page":"Locomotion API","title":"BehaviorDataNIR.get_reversal_events","text":"get_reversal_events(param, velocity, t_range, max_t)\n\nFinds reversal events.\n\nArguments:\n\nparam: Dictionary containing the following variables:\nrev_len_thresh: Number of consecutive reversal time points necessary for a reversal event\nrev_v_thresh: Velocity threshold below which the worm is counted as reversing\nvelocity: Worm velocity\nt_range: Time range over which to compute reversal events\nmax_t: Maximum time point\n\n\n\n\n\n","category":"function"},{"location":"locomotion/#BehaviorDataNIR.compute_reversal_times","page":"Locomotion API","title":"BehaviorDataNIR.compute_reversal_times","text":"compute_reversal_times(reversals, max_t)\n\nComputes the duration of each reversal event.\n\nArguments:\n\nreversals: List of all time points where the worm is reversing\nmax_t: Maximum time point in dataset\n\n\n\n\n\n","category":"function"},{"location":"locomotion/#Stage-data","page":"Locomotion API","title":"Stage data","text":"","category":"section"},{"location":"locomotion/","page":"Locomotion API","title":"Locomotion API","text":"zero_stage","category":"page"},{"location":"locomotion/#BehaviorDataNIR.zero_stage","page":"Locomotion API","title":"BehaviorDataNIR.zero_stage","text":"zero_stage(x::Array{<:AbstractFloat,2})\n\nSets the first time point to start at 0, 0\n\nArguments\n\npos_stage: (x,y) location, 2 by T array where T is len(time points)\n\n\n\n\n\n","category":"function"},{"location":"#BehaviorDataNIR.jl-Documentation","page":"BehaviorDataNIR.jl Documentation","title":"BehaviorDataNIR.jl Documentation","text":"","category":"section"},{"location":"","page":"BehaviorDataNIR.jl Documentation","title":"BehaviorDataNIR.jl Documentation","text":"The BehaviorDataNIR.jl package provides a collection of utilities for computing behavior information from NIR recordings.","category":"page"},{"location":"","page":"BehaviorDataNIR.jl Documentation","title":"BehaviorDataNIR.jl Documentation","text":"Pages = [\"locomotion.md\", \"posture.md\", \"feeding.md\", \"data.md\", \"sync.md\", \"util.md\"]","category":"page"},{"location":"util/#Utilities-API","page":"Utilities API","title":"Utilities API","text":"","category":"section"},{"location":"util/#Angle-utilities","page":"Utilities API","title":"Angle utilities","text":"","category":"section"},{"location":"util/","page":"Utilities API","title":"Utilities API","text":"recenter_angle\nlocal_recenter_angle\nvec_to_angle\nmake_vec\nang_btw_vec","category":"page"},{"location":"util/#BehaviorDataNIR.recenter_angle","page":"Utilities API","title":"BehaviorDataNIR.recenter_angle","text":"recenter_angle(angle; ref=0)\n\nRecenters angle to be within pi of a reference angle ref (optional, default 0)\n\n\n\n\n\n","category":"function"},{"location":"util/#BehaviorDataNIR.local_recenter_angle","page":"Utilities API","title":"BehaviorDataNIR.local_recenter_angle","text":"local_recenter_angle(angles; delta=10)\n\nRecenters angles to be continuous. delta (default 10) is the timespan of reference angles.\n\n\n\n\n\n","category":"function"},{"location":"util/#BehaviorDataNIR.vec_to_angle","page":"Utilities API","title":"BehaviorDataNIR.vec_to_angle","text":"vec_to_angle(vec)\n\nConverts a vector into an angle in the lab xy-coordinate space.\n\n\n\n\n\n","category":"function"},{"location":"util/#BehaviorDataNIR.make_vec","page":"Utilities API","title":"BehaviorDataNIR.make_vec","text":"make_vec(x::Array{<:AbstractFloat,1}, y::Array{<:AbstractFloat,1})\n\nCreates a vector out of x and y position variables.\n\n\n\n\n\n","category":"function"},{"location":"util/#Other-utilities","page":"Utilities API","title":"Other utilities","text":"","category":"section"},{"location":"util/","page":"Utilities API","title":"Utilities API","text":"impute_list\nget_lsqerr\nsavitzky_golay_filter\neuclidean_dist","category":"page"},{"location":"util/#BehaviorDataNIR.impute_list","page":"Utilities API","title":"BehaviorDataNIR.impute_list","text":"impute_list(x::Array{<:AbstractFloat,1})\n\nImputes missing data (NaN or missing) with interpolation\n\nArguments\n\nx: 1D data to impute\n\n\n\n\n\n","category":"function"},{"location":"util/#BehaviorDataNIR.get_lsqerr","page":"Utilities API","title":"BehaviorDataNIR.get_lsqerr","text":"get_lsqerr(fit, raw)\n\nComputes least squares error of a fit.\n\n\n\n\n\n","category":"function"},{"location":"util/#BehaviorDataNIR.savitzky_golay_filter","page":"Utilities API","title":"BehaviorDataNIR.savitzky_golay_filter","text":"savitzky_golay_filter(data, lag; is_derivative::Bool=false, has_inflection::Bool=true, is_angle=false)\n\nFilters data using the Savitzky-Golay algorithm, either for smoothing or differentiating purposes.\n\nArguments:\n\ndata: Data to filter.\nlag: Interval on each side to use for filtering.\nis_derivative::Bool (optional, default false): Whether to differentiate the data (vs just smoothing)\nhas_inflection::Bool (optional, default true): Whether the data has inflection points within smoothing interval.   If set to true, apply higher-order smoothing function.\nis_angle::Bool (optional, default false): Whether the data is an angle (ie -pi=pi)\n\n\n\n\n\n","category":"function"},{"location":"posture/#Worm-Posture-API","page":"Worm Posture API","title":"Worm Posture API","text":"","category":"section"},{"location":"posture/#Spline-computation","page":"Worm Posture API","title":"Spline computation","text":"","category":"section"},{"location":"posture/","page":"Worm Posture API","title":"Worm Posture API","text":"segment_worm!\ncompute_worm_spline!\ncompute_worm_thickness\nget_segment_end_matrix\ninterpolate_splines!","category":"page"},{"location":"posture/#BehaviorDataNIR.segment_worm!","page":"Worm Posture API","title":"BehaviorDataNIR.segment_worm!","text":"segment_worm!(model, img_raw, img_label; θ=0.75)\n\nSegment the worm given img (size: (968, 732)).\n\nArguments\n\nmodel: unet model\nimg_raw: raw img\nimg_label: pre-allocated label array <:Int32. Size: (480,360)\n\nOptional arguments\n\nθ=0.75: unet output threshold\n\n\n\n\n\n","category":"function"},{"location":"posture/#BehaviorDataNIR.compute_worm_spline!","page":"Worm Posture API","title":"BehaviorDataNIR.compute_worm_spline!","text":"compute_worm_spline!(\n    param, path_h5, worm_seg_model, worm_thickness, med_axis_dict, pts_order_dict, is_omega_dict,\n    x_array, y_array, nir_worm_angle, eccentricity; timepts=\"all\"\n)\n\nCompute the worm spline for a given set of parameters. Writes to most of its input parameters.\n\nArguments:\n\nparam: A dictionary containing the parameters for the worm spline computation.\npath_h5: The path to the HDF5 file containing the worm data.\nworm_seg_model: The worm segmentation model.\nworm_thickness: The thickness of the worm.\nmed_axis_dict: A dictionary containing where to write medial axis data.\npts_order_dict: A dictionary containing where to write the order of spline points.\nis_omega_dict: A dictionary containing where to write values of whether the worm is self-intersecting.\nx_array: An array containing the x-coordinates of the worm splines. Will be modified.\ny_array: An array containing the y-coordinates of the worm splines. Will be modified.\nnir_worm_angle: The angle of the worm.\neccentricity: The eccentricity of the worm.\ntimepts: The timepoints to compute the worm spline for. Defaults to \"all\".\n\nReturns:\n\nimg_label: An array containing the labeled image.\nerrors: A dictionary containing any errors that occurred during computation.\n\n\n\n\n\n","category":"function"},{"location":"posture/#BehaviorDataNIR.compute_worm_thickness","page":"Worm Posture API","title":"BehaviorDataNIR.compute_worm_thickness","text":"compute_worm_thickness(param, path_h5, worm_seg_model, med_axis_dict, is_omega_dict)\n\nThis function computes the thickness of the worm at each point along its length.\n\nArguments:\n\nparam: A dictionary containing parameters for the function.\npath_h5: The path to the HDF5 file containing the worm images.\nworm_seg_model: The segmentation model used to segment the worm images.\nmed_axis_dict: A dictionary containing the medial axis of the worm at each timepoint.\nis_omega_dict: A dictionary indicating whether each timepoint corresponds to an omega turn.\n\nReturns:\n\ndists: A vector containing the thickness of the worm at each point along its length.\n\n\n\n\n\n","category":"function"},{"location":"posture/#BehaviorDataNIR.get_segment_end_matrix","page":"Worm Posture API","title":"BehaviorDataNIR.get_segment_end_matrix","text":"get_segment_end_matrix(param, x_array, y_array)\n\nComputes equally-spaced points along the worm spline.\n\nArguments:\n\nparam: Dictionary of parameters that includes the following values:\nnum_center_pts: Number of points along the spline\nsegment_len: Length of each segment\nx_array: x-locations of spline at each time point\ny_array: y-locations of spline at each time point\n\n\n\n\n\n","category":"function"},{"location":"posture/#BehaviorDataNIR.interpolate_splines!","page":"Worm Posture API","title":"BehaviorDataNIR.interpolate_splines!","text":"interpolatesplines!(datadict)\n\nInterpolates worm splines to time points where the spline computation crashed.\n\n\n\n\n\n","category":"function"},{"location":"posture/#Curvature-computation","page":"Worm Posture API","title":"Curvature computation","text":"","category":"section"},{"location":"posture/","page":"Worm Posture API","title":"Worm Posture API","text":"get_body_angles!\nget_angular_velocity!\nget_curvature_variables!\nget_nose_curling!\nangular_velocity\nget_worm_body_angle\nget_worm_vector\nget_tot_worm_curvature\nnose_curling","category":"page"},{"location":"posture/#BehaviorDataNIR.get_body_angles!","page":"Worm Posture API","title":"BehaviorDataNIR.get_body_angles!","text":"getbodyangles!(data_dict::Dict, param::Dict; prefix::String=\"\")\n\nAdds body angles to data_dict given param. Can add a prefix (default empty string) to all confocal variables.\n\n\n\n\n\n","category":"function"},{"location":"posture/#BehaviorDataNIR.get_angular_velocity!","page":"Worm Posture API","title":"BehaviorDataNIR.get_angular_velocity!","text":"get_angular_velocity!(data_dict::Dict, param::Dict; prefix::String=\"\")\n\nGets angular velocity from data_dict and param. Can add a prefix (default empty string) to all confocal variables.\n\n\n\n\n\n","category":"function"},{"location":"posture/#BehaviorDataNIR.get_curvature_variables!","page":"Worm Posture API","title":"BehaviorDataNIR.get_curvature_variables!","text":"get_curvature_variables!(data_dict::Dict, param::Dict; prefix::String=\"\")\n\nGets curvature, head angle, and related variables from data_dict and param. Can add a prefix (default empty string) to all confocal variables.\n\n\n\n\n\n","category":"function"},{"location":"posture/#BehaviorDataNIR.get_nose_curling!","page":"Worm Posture API","title":"BehaviorDataNIR.get_nose_curling!","text":"get_nose_curling!(data_dict::Dict, param::Dict; prefix::String=\"\")\n\nGets self intersection variables from data_dict and param. Can add a prefix (default empty string) to all confocal variables.\n\n\n\n\n\n","category":"function"},{"location":"posture/#BehaviorDataNIR.get_worm_body_angle","page":"Worm Posture API","title":"BehaviorDataNIR.get_worm_body_angle","text":"get_worm_body_angle(x_array, y_array, segment_end_matrix, pts)\n\nGets the body angle of a worm between the given points in the spline.\n\nArguments:\n\nx_array: Segmentation array in the x-dimension\ny_array: Segmentation array in the y-dimension\nsegment_end_matrix: Matrix of consistently-spaced segmentation locations across time points\npts: A list of three points in the spline to get the body angle between.\n\n\n\n\n\n","category":"function"},{"location":"posture/#BehaviorDataNIR.get_worm_vector","page":"Worm Posture API","title":"BehaviorDataNIR.get_worm_vector","text":"get_worm_vector(x_array, y_array, segment_end_matrix, seg_range)\n\nGets the vector the worm is facing, in the lab xy-coordinate system.\n\nArguments:\n\nx_array: Segmentation array in the x-dimension\ny_array: Segmentation array in the y-dimension\nsegment_end_matrix: Matrix of consistently-spaced segmentation locations across time points\nseg_range: Set of segmentation locations to use to compute the centroid (which will be the average of them)\n\n\n\n\n\n","category":"function"},{"location":"posture/#BehaviorDataNIR.get_tot_worm_curvature","page":"Worm Posture API","title":"BehaviorDataNIR.get_tot_worm_curvature","text":"get_tot_worm_curvature(body_angle, min_len; directional::Bool=false)\n\nComputes total worm curvature\n\nArguments:\n\nbody_angle: Array of worm body angles\nmin_len: If there aren't this many angles at a given time point, interpolate that time point instead of computing it\ndirectional::Bool (default false): Use directional curvature.\n\n\n\n\n\n","category":"function"},{"location":"posture/#BehaviorDataNIR.nose_curling","page":"Worm Posture API","title":"BehaviorDataNIR.nose_curling","text":"nose_curling(spline_x, spline_y, segment_end; max_i=1)\n\nGets the smallest ratio between distance between two points in space and distance between them along the worm's curvature. A value of 0 means that the worm is intersecting itself, while a value of 1 means the worm is a straight line.\n\nArguments:\n\nspline_x: x positions in worm spline\nspline_y: y positions in worm spline\nsegment_end: equally-spaced segments along worm spline\nmax_i (default 1): Maximum location along the medium axis to try. For nose curling, use 1.\n\n\n\n\n\n","category":"function"}]
}
