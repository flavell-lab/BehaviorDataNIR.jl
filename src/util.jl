function diff_lag(A::AbstractVector, lag::Int)
    A[lag+1:end] .- A[1:end-lag]
end

function vec_ang(v1, v2)
    acos(clamp(dot(v1, v2) / (norm(v1, 2) * norm(v2, 2)), -1, 1))
end

function read_pos_feature(h5f::HDF5.File)
    pos_feature = read(h5f, "pos_feature")
    pos_feature_unet = copy(pos_feature)
    pos_feature_unet[:,1:2,:] = round.((pos_feature_unet[:,1:2,:] .- 4) ./ 2)
   
    pos_feature, pos_feature_unet
end

function read_pos_feature(path_h5)
    h5open(path_h5, "r") do h5f
        import_pos_feature(h5f)
    end
end

function read_h5(path_h5::String)
    h5open(path_h5, "r") do h5f
        read(h5f, "img_nir"), read(h5f, "pos_stage"), read_pos_feature(h5f)
    end
end