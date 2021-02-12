function nmp_vec(pos_feature)
    n = pos_feature[1,1:2,:] # nose
    m = pos_feature[2,1:2,:] # metacorpus
    p = pos_feature[3,1:2,:] # pharynx
    
    mn = n .- m
    mp = p .- m
    mp⊥ = Array(hcat(mp[2,:], - mp[1,:])')
    
    mn, mp, mp⊥
end