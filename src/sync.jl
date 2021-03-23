function detect_nir_timing(di_nir)
        # behavior camera - FLIR
    list_nir_on = findall(diff(di_nir) .> 1) .+ 1
    list_nir_off = findall(diff(di_nir) .< -1) .+ 1
    nir_record_on = diff(list_nir_on) .> 500
    nir_record_off = diff(list_nir_off) .> 500

    if list_nir_on[1] > 500 # no trigger before the first
        s_nir_start = list_nir_on[1]
    elseif sum(nir_record_on) == 2
        s_nir_start = list_nir_on[findfirst(nir_record_on) + 1]
    else
        error("more than 2 recording on detected for FLIR camera")
    end

    if list_nir_off[end] < length(di_nir) - 500
        s_nir_stop = list_nir_off[end]
    elseif sum(nir_record_off) <= 2
        s_nir_stop = list_nir_off[findlast(diff(list_nir_off) .> 500)] 
    else
        error("more than 2 recording off detected for FLIR camera")
    end

    list_nir_on = filter(x -> s_nir_start - 5 .< x .< s_nir_stop .+ 5, list_nir_on)
    list_nir_off = filter(x -> s_nir_start - 5 .< x .< s_nir_stop .+ 5, list_nir_off)
    
    if length(list_nir_on) != length(list_nir_off)
        error("length(list_nir_on) != length(list_nir_off)")
    end
    
    list_nir_on, list_nir_off
end

function detect_confocal_timing(ai_laser)
    ai_laser_bin = Int16.(ai_laser .> mean(ai_laser)) # binarize laser analog signal

    list_confocal_on = findall(diff(ai_laser_bin) .== 1) .+ 1
    list_confocal_off = findall(diff(ai_laser_bin) .== -1) .+ 1

    list_stack_start = list_confocal_on[findall(diff(list_confocal_on) .> 300) .+ 1]
    prepend!(list_stack_start, list_confocal_on[1])
    list_stack_stop = list_confocal_off[findall(diff(list_confocal_off) .> 300)]
    append!(list_stack_stop, list_confocal_off[end])

    if length(list_stack_start) != length(list_stack_stop)
        error("n(stack_off_confocal) != n(stack_on_confocal)")
    end
    
    list_stack_start, list_stack_stop
end

function sync_timing(di_nir, ai_laser, img_id, q_iter_save, n_img_nir)
    list_nir_on, list_nir_off = detect_nir_timing(di_nir)
    list_stack_start, list_stack_stop = detect_confocal_timing(ai_laser)
    
    img_id_diff = diff(img_id)
    prepend!(img_id_diff, 1)
    if length(list_nir_on) - sum(diff(img_id)) > 3
        error("the detected trigger count is different from the image id data by more than 3")
    else
        img_id_diff[end] += length(list_nir_on) - sum(diff(img_id)) - 1
    end

    idx_nir_save = Bool[]
    for (Δn, q_save) = zip(img_id_diff, q_iter_save)
        if Δn == 1
            push!(idx_nir_save, q_save)
        else
            push!(idx_nir_save, q_save)
            for i = 1:Δn-1
                push!(idx_nir_save, false)
            end
        end
    end

    if sum(idx_nir_save) != n_img_nir
        error("detected number of NIR frames != saved NIR frames")
    end

    timing_stack = hcat(list_stack_start, list_stack_stop)
    timing_nir = hcat(list_nir_on, list_nir_off)[idx_nir_save,:]

    confocal_to_nir = []
    nir_to_confocal = zeros(size(timing_nir,1))
    for i = 1:size(timing_stack, 1)
        start_, end_ = timing_stack[i,:]

        nir_on_bit = start_ .< timing_nir[:,1] .< end_
        nir_off_bit = start_ .< timing_nir[:,2] .< end_

        idx_ = findall(nir_on_bit .& nir_off_bit)
        push!(confocal_to_nir, idx_)

        nir_to_confocal[idx_] .= i
    end
        
    confocal_to_nir, nir_to_confocal, timing_stack, timing_nir
end

function sync_timing(path_h5)
    n_img_nir, daqmx_ai, daqmx_di, img_metadata = h5open(path_h5, "r") do h5f
        n_img_nir = size(h5f["img_nir"])[3]
        daqmx_ai = read(h5f, "daqmx_ai")
        daqmx_di = read(h5f, "daqmx_di")
        img_metadata = read(h5f, "img_metadata")
        n_img_nir, daqmx_ai, daqmx_di, img_metadata
    end

    ai_laser = daqmx_ai[:,1]
    ai_piezo = daqmx_ai[:,2]
    ai_stim = daqmx_ai[:,3]
    di_confocal = Float32.(daqmx_di[:,1])
    di_nir = Float32.(daqmx_di[:,2])
    img_timestamp = img_metadata["img_timestamp"]
    img_id = img_metadata["img_id"]
    q_iter_save = img_metadata["q_iter_save"]
    
    sync_timing(di_nir, ai_laser, img_id, q_iter_save, n_img_nir)
end

function sync_stim(stim, timing_stack, timing_nir)
    stim_to_confocal = [mean(stim[timing_stack[i,1]:timing_stack[i,2]]) for i = 1:size(timing_stack,1)]
    stim_to_nir = stim[round.(Int, dropdims(mean(timing_nir, dims=2), dims=2))]

    stim_to_confocal, stim_to_nir
end

function signal_stack_repeatability(signal, timing_stack; sampling_rate=5000)
    s_stack_start = timing_stack[:,1]
    s_stack_end = timing_stack[:,2]
    n_stack = size(s_stack_end, 1)
    n_stack_len = minimum(s_stack_end - s_stack_start)
    signal_eta = zeros(n_stack, n_stack_len)
    for i = 1:n_stack
        signal_eta[i,:] .= signal[s_stack_start[i]:s_stack_start[i]+n_stack_len-1]
    end

    signal_eta_u = dropdims(mean(signal_eta, dims=1), dims=1)
    signal_eta_s = dropdims(std(signal_eta, dims=1), dims=1)
    list_t = collect(1:n_stack_len) / sampling_rate
    
    signal_eta_u, signal_eta_s, list_t, n_stack
end