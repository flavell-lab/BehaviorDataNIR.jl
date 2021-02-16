function zero_stage(pos_stage)
    x, y = pos_stage[1,:], pos_stage[2,:]
    @assert(all(.!isnan.([x[1], y[1]])))

    x .- x[1], y .- y[1]
end

"""
    impute_stage(x::Array{<:AbstractFloat,1})

Sets the first time point to start at 0, 0

Arguments
---------
* `pos_stage`: (x,y) location, 2 by T array where T is len(time points)
"""
function impute_stage(x::Array{<:AbstractFloat,1})
    x, y = pos_stage[1,:], pos_stage[2,:]
    @assert(all(.!isnan.([x[1], y[1]])))
    
    x .- x[1], y .- y[1]
end

"""
    impute_stage(x::Array{<:AbstractFloat,1})

Imputes missing data (NaN or missing) with interpolation

Arguments
---------
* `x`: 1D data to impute
"""
function impute_stage(x::Array{<:AbstractFloat,1})
    convert.(eltype(x), Impute.interp(replace(x, NaN=>missing)))
end

"""
    speed(Δx::Array{<:AbstractFloat,1}, Δy::Array{<:AbstractFloat,1}, Δt::AbstractFloat)

Computes speed

Arguments
---------
* `Δx`: discrete difference of x
* `Δy`: discrete difference of y
* `Δt`: time interval
"""
function speed(Δx::Array{<:AbstractFloat,1}, Δy::Array{<:AbstractFloat,1}, Δt::AbstractFloat)
    sqrt.(Δx .^2 .+ Δy .^2) / Δt
end

"""
    speed(x, y; lag::Int, fps=FLIR_FPS)

Computes speed using x, y coordinates

`lag` determines the interval at which to compute the discrete difference

Arguments
---------
* `x`: list of x
* `y`: list of y
* `lag`: number of time points for discerete difference interval
* `fps`: fps
"""
function speed(x, y; lag::Int, fps=FLIR_FPS)
    Δx, Δy = diff_lag.([x, y], lag=lag)
    Δt = lag * 1 / fps
    
    speed(Δx, Δy, Δt)
end

function Δpos_angle(Δx, Δy)
    atan.(Δy ./ Δx)
end

function Δpos_angle(x, y; lag::Int)
    Δx, Δy = diff_lag.([x, y], lag=lag)
    
    Δpos_angle(Δx, Δy)
end

function angular_velocity(Δθ, Δt)
    Δθ / Δt # rad/s
end

function angular_velocity(x, y; lag::Int, fps=FLIR_FPS)
    Δx, Δy = diff_lag.([x, y], lag=lag)
    Δθ = diff(Δpos_angle(Δx, Δy))
    Δt = (1 / fps) * lag
    
    angular_velocity(Δθ, Δt) # rad/s
end
