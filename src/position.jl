function zero_stage(pos_stage)
    x, y = pos_stage[1,:], pos_stage[2,:]
    @assert(all(.!isnan.([x[1], y[1]])))
    
    x .- x[1], y .- y[1]
end

function impute_stage(x)
    convert.(eltype(x), Impute.interp(replace(x, NaN=>missing)))
end