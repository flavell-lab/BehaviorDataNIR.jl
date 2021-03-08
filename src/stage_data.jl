"""
    zero_stage(x::Array{<:AbstractFloat,2})

Sets the first time point to start at 0, 0

Arguments
---------
* `pos_stage`: (x,y) location, 2 by T array where T is len(time points)
"""
function zero_stage(pos_stage::Array{<:AbstractFloat,2})
    x, y = pos_stage[1,:], pos_stage[2,:]
    @assert(all(.!isnan.([x[1], y[1]])))

    x .- x[1], y .- y[1]
end


"""
    impute_list(x::Array{<:AbstractFloat,1})

Imputes missing data (NaN or missing) with interpolation

Arguments
---------
* `x`: 1D data to impute
"""
function impute_list(x::Array{<:AbstractFloat,1})
    imputed_lst = Impute.interp(replace(x, NaN=>missing))
    if typeof(imputed_lst[1]) == Missing
        imputed_lst = impute(imputed_lst, Impute.NOCB())
    end
        
    if typeof(imputed_lst[end]) == Missing
        imputed_lst = impute(imputed_lst, Impute.LOCF())
    end
    
    convert.(eltype(x), imputed_lst)
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
function speed(x, y, lag::Int; fps=FLIR_FPS)
    Δx, Δy = diff_lag.([x, y], lag)
    Δt = lag * 1 / fps
    
    speed(Δx, Δy, Δt)
end

function time_axis(list::AbstractVector, lag=0; fps=FLIR_FPS)
    num_frame = maximum(size(list))
    collect(1 : num_frame - lag) * 1 / fps
end

function Δpos_angle(Δx, Δy)
    atan.(Δy ./ Δx)
end

function Δpos_angle(x, y, lag::Int)
    Δx, Δy = diff_lag.([x, y], lag) # diff_lag in util.jl
    
    Δpos_angle(Δx, Δy)
end

function angular_velocity(Δθ, Δt)
    Δθ / Δt # rad/s
end

function angular_velocity(x, y, lag::Int; fps=FLIR_FPS)
    Δx, Δy = diff_lag.([x, y], lag)
    Δθ = diff(Δpos_angle(Δx, Δy))
    Δt = (1 / fps) * lag
    
    angular_velocity(Δθ, Δt) # rad/s
end

function ang_btw_vec(v1::Array{<:AbstractFloat,2}, v2::Array{<:AbstractFloat,2})
    timept = min(size(v1, 2), size(v2, 2))
    vec_ang_lst = zeros(timept)
    
    for i = 1 : timept
        try
            vec_ang_lst[i] = vec_ang(v1[:, i], v2[:, i])
        catch
            @warn "error at index " * string(i) * "// v1: " * string(v1[:, i]) * "// v2: " * string(v2[:, i])
        end
    end
    
    vec_ang_lst
end

# take in 2 x num_vec array as iput and return 1 x num_vec array of magnitudes
function magnitude_vec(list_v::Array{<:AbstractFloat,2})
    sqrt.(list_v[1, :] .^ 2 .+ list_v[2, :] .^ 2)
end

# A = [0.9 0.1; 0.9 0.1] # transition matrix
# B = [Normal(170,10), Normal(13,5)] # observation dist, mean and variance
function reversal_state(list::Array{<:AbstractFloat,1}, A, B)
    hmm = HMM(A, B)
    viterbi(hmm, list)
end

function make_vec(x::Array{<:AbstractFloat,1}, y::Array{<:AbstractFloat,1})
    convert(Array{Float64,2}, (hcat(x, y))') # can't do typeof(x) since it's 1D
end