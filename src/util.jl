function diff_lag(A::AbstractVector; lag::Int)
    A[lag+1:end] .- A[1:end-lag]
end

function vec_ang(v1, v2)
    acos(dot(v1, v2) / (norm(v1,2) * norm(v2,2)))
end

function read_hdf5(input_h5_path::String)
    img_nir = h5read(input_h5_path, "img_nir")
    list_pos_stage = h5read(input_h5_path, "pos_stage")
    pos_feature = h5read(input_h5_path, "pos_feature")
    
    return img_nir, list_pos_stage, pos_feature
end