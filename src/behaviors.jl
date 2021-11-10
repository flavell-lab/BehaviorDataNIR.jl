"""
Adds body angles to `data_dict` given `param`.
"""
function get_body_angles!(data_dict::Dict, param::Dict)
    s = size(data_dict["x_array"],1)
    m = maximum([length(x) for x in data_dict["segment_end_matrix"]])-1
    data_dict["nir_body_angle"] = zeros(param["max_pt"]-1,s)
    data_dict["nir_body_angle_all"] = zeros(m,s)
    data_dict["nir_body_angle_absolute"] = zeros(m,s)
    data_dict["body_angle"] = zeros(param["max_pt"]-1,maximum(param["t_range"]))
    data_dict["body_angle_all"] = zeros(m,maximum(param["t_range"]))
    data_dict["body_angle_absolute"] = zeros(m,maximum(param["t_range"]))
    for pos in 1:m
        for t in 1:s
            if length(data_dict["segment_end_matrix"][t]) > pos
                Δx = data_dict["x_array"][t,data_dict["segment_end_matrix"][t][pos+1]] - data_dict["x_array"][t,data_dict["segment_end_matrix"][t][pos]]
                Δy = data_dict["y_array"][t,data_dict["segment_end_matrix"][t][pos+1]] - data_dict["y_array"][t,data_dict["segment_end_matrix"][t][pos]]
                data_dict["nir_body_angle_absolute"][pos,t] = recenter_angle(vec_to_angle([Δx, Δy])[1])
                data_dict["nir_body_angle_all"][pos,t] = recenter_angle(vec_to_angle([Δx, Δy])[1] - data_dict["nir_worm_angle"][t])
            else
                data_dict["nir_body_angle_absolute"][pos,t]= NaN
                data_dict["nir_body_angle_all"][pos,t]= NaN
            end
        end
        data_dict["nir_body_angle_absolute"][pos,:] .= local_recenter_angle(data_dict["nir_body_angle_absolute"][pos,:], delta=param["body_angle_t_lag"])
    end
    
    for t in 1:s
        data_dict["nir_body_angle_absolute"][:,t] .= local_recenter_angle(data_dict["nir_body_angle_absolute"][:,t], delta=param["body_angle_pos_lag"])
        data_dict["nir_body_angle_all"][:,t] .= local_recenter_angle(data_dict["nir_body_angle_all"][:,t], delta=param["body_angle_pos_lag"])
    end
    
    for pos in 1:m
        data_dict["body_angle_all"][pos,:] .= vec_to_confocal(data_dict["nir_body_angle_all"][pos,:])
        data_dict["body_angle_absolute"][pos,:] .= vec_to_confocal(data_dict["nir_body_angle_absolute"][pos,:])
        
        if pos < param["max_pt"]
            data_dict["nir_body_angle"][pos,:] .= data_dict["nir_body_angle_all"][pos,:]      
            data_dict["nir_body_angle"][pos,:] .= impute_list(data_dict["nir_body_angle"][pos,:])
            data_dict["body_angle"][pos,:] .= vec_to_confocal(data_dict["nir_body_angle"][pos,:])
            data_dict["body_angle_absolute"][pos,:] .= vec_to_confocal(data_dict["nir_body_angle_absolute"][pos,:]);
        end
    end
end

"""
Gets angular velocity from `data_dict` and `param`.
"""
function get_angular_velocity!(data_dict::Dict, param::Dict)
    data_dict["worm_angle"] = vec_to_confocal(data_dict["nir_worm_angle"])

    turning_angle = impute_list(data_dict["body_angle_absolute"][param["head_pts"][1],:])

    data_dict["worm_head_facing_angle_filtered"] = gstv(turning_angle .- turning_angle[1], m, l);
    data_dict["angular_velocity"] = zeros(param["max_t"])
    data_dict["angular_velocity"][2:end-1] = diff_lag(data_dict["worm_head_facing_angle_filtered"], 2)
    data_dict["angular_velocity"][1] = data_dict["angular_velocity"][2]
    data_dict["angular_velocity"][end-1] = data_dict["angular_velocity"][end];
end

"""
Gets velocity, speed, reversal, and related variables from `data_dict` and `param`.
"""
function get_velocity!(data_dict::Dict, param::Dict)
    data_dict["filt_xmed"] = gstv(Float64.(data_dict["x_med"]), param["v_stage_m_filt"], param["v_stage_λ_filt"])
    data_dict["filt_ymed"] = gstv(Float64.(data_dict["y_med"]), param["v_stage_m_filt"], param["v_stage_λ_filt"]);

    Δx = [0.0]
    Δy = [0.0]
    append!(Δx, diff(data_dict["filt_xmed"]))
    append!(Δy, diff(data_dict["filt_ymed"]))
    Δt = 1.0 / param["FLIR_FPS"]
    data_dict["nir_mov_vec_stage"] = make_vec(Δx, Δy)
    data_dict["mov_vec_stage"] = vec_to_confocal(data_dict["nir_mov_vec_stage"])
    data_dict["nir_mov_angle_stage"] = impute_list(vec_to_angle(data_dict["nir_mov_vec_stage"]))
    data_dict["mov_angle_stage"] = vec_to_confocal(data_dict["nir_mov_angle_stage"])
    data_dict["nir_speed_stage"] = speed(Δx, Δy, Δt)
    data_dict["speed_stage"] = vec_to_confocal(data_dict["nir_speed_stage"])
    mn_vec, mp_vec, orthog_mp_vec = nmp_vec(pos_feature);
    data_dict["pm_angle"] = vec_to_angle(mp_vec);
    data_dict["nir_velocity_stage"] = data_dict["nir_speed_stage"] .* cos.(data_dict["nir_mov_angle_stage"] .- data_dict["pm_angle"])
    data_dict["velocity_stage"] = vec_to_confocal(data_dict["nir_velocity_stage"])
    data_dict["reversal_events"], data_dict["all_rev"] = get_reversal_events(param, data_dict["velocity_stage"], param["t_range"]);
    data_dict["rev_times"] = compute_reversal_times(data_dict["all_rev"], maximum(param["t_range"]));
end

"""
Gets curvature, head angle, and related variables from `data_dict` and `param`.
"""
function get_curvature_variables!(data_dict::Dict, param::Dict)
    data_dict["worm_curvature"] = get_tot_worm_curvature(data_dict["body_angle"], size(data_dict["body_angle"],1));
    data_dict["ventral_worm_curvature"] = get_tot_worm_curvature(data_dict["body_angle"], size(data_dict["body_angle"],1), directional=true);
    data_dict["nir_head_angle"] = -get_worm_body_angle(data_dict["x_array"], data_dict["y_array"], data_dict["segment_end_matrix"], param["head_pts"])
    data_dict["nir_nose_angle"] = -get_worm_body_angle(data_dict["x_array"], data_dict["y_array"], data_dict["segment_end_matrix"], param["nose_pts"])

    data_dict["head_angle"] = vec_to_confocal(data_dict["nir_head_angle"])
    data_dict["nose_angle"] = vec_to_confocal(data_dict["nir_nose_angle"]);
end

"""
Gets self intersection variables from `data_dict` and `param`.
"""
function get_self_intersection!(data_dict::Dict, param::Dict)
    data_dict["nir_self_intersect_ratio"] = Vector{Float64}()
    data_dict["nir_self_intersect_ratio_head"] = Vector{Float64}()
    max_med_len = param["segment_len"] * param["max_pt"]
    for t=1:data_dict["max_t_nir"]
        if length(data_dict["med_axis_dict"][t][1]) < max_med_len
            push!(data_dict["nir_self_intersect_ratio"], NaN)
            push!(data_dict["nir_self_intersect_ratio_head"], NaN)
        else
            push!(data_dict["nir_self_intersect_ratio"], 1 ./self_intersect_ratio(map(x->x[1:max_med_len], data_dict["med_axis_dict"][t])))
            push!(data_dict["nir_self_intersect_ratio_head"], 1 ./self_intersect_ratio(map(x->x[1:max_med_len], data_dict["med_axis_dict"][t]), max_i=param["segment_len"]*param["head_pts"][3]))
        end
    end
    data_dict["nir_self_intersect_ratio"] = impute_list(data_dict["nir_self_intersect_ratio"])
    data_dict["nir_self_intersect_ratio_head"] = impute_list(data_dict["nir_self_intersect_ratio_head"]);
    data_dict["self_intersect_ratio"] = vec_to_confocal(data_dict["nir_self_intersect_ratio"])
    data_dict["self_intersect_ratio_head"] = vec_to_confocal(data_dict["nir_self_intersect_ratio_head"]);
end

"""
Modifies `combined_data_dict` to add merged NIR data from two datasets `data_dict` and `data_dict_2`,
together with a parameter file `param` for the primary dataset.
"""
function merge_nir_data!(combined_data_dict::Dict, data_dict::Dict, data_dict_2::Dict, param::Dict)
    combined_data_dict["confocal_to_nir_1"] = data_dict["confocal_to_nir"]
    combined_data_dict["confocal_to_nir_2"] = data_dict_2["confocal_to_nir"]
    combined_data_dict["nir_to_confocal_1"] = data_dict["nir_to_confocal"]
    combined_data_dict["nir_to_confocal_2"] = data_dict_2["nir_to_confocal"]
    combined_data_dict["max_t_nir_1"] = data_dict["max_t_nir"]
    combined_data_dict["max_t_nir_2"] = data_dict_2["max_t_nir"]


    for var in param["concat_vars"]
        if length(size(data_dict[var])) == 1
            combined_data_dict[var] = zeros(data_dict["max_t_all"])
            combined_data_dict[var][1:length(data_dict[var])] .= data_dict[var]
            combined_data_dict[var][length(data_dict[var])+1:end] .= data_dict_2[var]
        elseif length(size(data_dict[var])) == 2
            max_size = max(size(data_dict[var],1), size(data_dict_2[var],1))
            combined_data_dict[var] = zeros(max_size, data_dict["max_t_all"])
            combined_data_dict[var][1:size(data_dict[var],1),1:size(data_dict[var],2)] .= data_dict[var][:,:]
            combined_data_dict[var][1:size(data_dict_2[var],1),size(data_dict[var],2)+1:end] .= data_dict_2[var][:,:]
        else
            throw(ErrorException("number of dimensions must be 1 or 2"))
        end
    end
    for var in param["t_concat_vars"]
        combined_data_dict[var] = deepcopy(data_dict[var])
        append!(combined_data_dict[var], param["max_t"] .+ data_dict_2[var])
    end
    for var in param["nir_concat_vars"]
        if length(size(data_dict[var])) == 1
            tot_len = length(data_dict[var]) + length(data_dict_2[var])
            combined_data_dict[var] = zeros(tot_len)
            combined_data_dict[var][1:length(data_dict[var])] .= data_dict[var]
            combined_data_dict[var][length(data_dict[var])+1:end] .= data_dict_2[var]
        elseif length(size(data_dict[var])) == 2
            max_size = max(size(data_dict[var],1), size(data_dict_2[var],1))
            tot_len = size(data_dict[var],2) + size(data_dict_2[var],2)
            combined_data_dict[var] = zeros(max_size, tot_len)
            combined_data_dict[var][1:size(data_dict[var],1),1:size(data_dict[var],2)] .= data_dict[var][:,:]
            combined_data_dict[var][1:size(data_dict_2[var],1),size(data_dict[var],2)+1:end] .= data_dict_2[var][:,:]
        else
            throw(ErrorException("number of dimensions must be 1 or 2"))
        end
    end
end

"""
Import pumping data into a combined dataset from a csv file.
"""
function import_pumping!(combined_data_dict::Dict, path_pumping)
    combined_data_dict["pumping_nir"] = Float64[]
    combined_data_dict["pumping_raw"] = Float64[]
    combined_data_dict["pumping"] = Float64[]

    for (d, file) in enumerate(path_pumping)
        pumping = readdlm(file, ',', Any, '\n')
        pumping_nir = param["FLIR_FPS"] .* [(t in floor.(Int64.(pumping[2:end,2])./50)) ? 1 : 0 for t in 1:combined_data_dict["max_t_nir_$d"]]
        append!(combined_data_dict["pumping_nir"], pumping_nir)
        pumping_raw = nir_vec_to_confocal(pumping_nir, combined_data_dict["confocal_to_nir_$d"], length(combined_data_dict["confocal_to_nir_$d"]))
        if length(path_pumping == 1) 
            append!(pumping_raw, nir_vec_to_confocal(pumping_nir, combined_data_dict["confocal_to_nir_2"], length(combined_data_dict["confocal_to_nir_2"])))
        end
        append!(combined_data_dict["pumping_raw"], pumping_raw)
        append!(combined_data_dict["pumping"], gstv(pumping_raw, 10, 0.2))
    end
end
