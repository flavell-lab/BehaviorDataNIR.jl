function read_hdf5(input_h5_path::String)
    img_nir = h5read(input_h5_path, "img_nir")
    pos_stage = h5read(input_h5_path, "pos_stage")
    pos_feature = h5read(input_h5_path, "pos_feature")
    
    return img_nir, pos_stage, pos_feature
end

function nmp_vec(pos_feature)
    n = pos_feature[1,1:2,:] # nose
    m = pos_feature[2,1:2,:] # metacorpus
    p = pos_feature[3,1:2,:] # pharynx
    
    mn = n .- m
    mp = p .- m
    mp⊥ = Array(hcat(mp[2,:], - mp[1,:])')
    
    mn, mp, mp⊥
end