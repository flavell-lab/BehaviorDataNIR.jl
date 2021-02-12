function zero_stage(pos_stage)
    x, y = pos_stage[1,:], pos_stage[2,:]
    @assert(all(.!isnan.([x[1], y[1]])))
    
    x .- x[1], y .- y[1]
end

function impute_stage(x)
    convert.(eltype(x), Impute.interp(replace(x, NaN=>missing)))
end

function speed(x, y; lag::Int, fps=FLIR_FPS)
    Δx, Δy = diff_lag.([x, y], lag=lag)
    Δt = lag * 1 / fps
    
    sqrt.(Δx .^2 .+ Δy .^2) / Δt
end

function Δpos_angle(Δx, Δy)
    atan.(Δy ./ Δx)
end

function angular_velocity(x, y; lag::Int, fps=FLIR_FPS)
    Δx, Δy = diff_lag.([x, y], lag=lag)
    
    Δpos_angle(Δx, Δy) / ((1 / fps) * lag) # rad/s
end
