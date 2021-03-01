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
    imputed_lst = impute(replace(x, NaN => missing), Impute.Interpolate(); dims=:rows))
    
    if typeof(imputed_lst[1]) == Missing
        imputed_lst = impute(imputed_lst, Impute.NOCB(); dims=:rows)
    end
        
    if typeof(imputed_lst[end]) == Missing
        imputed_lst = impute(imputed_lst, Impute.LOCF(); dims=:rows)
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
function speed(x, y; lag::Int, fps=FLIR_FPS)
    Δx, Δy = diff_lag.([x, y], lag=lag)
    Δt = lag * 1 / fps
    
    speed(Δx, Δy, Δt)
end

function time_axis(list::AbstractVector; lag=0, fps=FLIR_FPS)
    num_frame = maximum(size(list))
    collect(1 : num_frame - lag) * 1 / fps
end

function Δpos_angle(Δx, Δy)
    atan.(Δy ./ Δx)
end

function Δpos_angle(x, y; lag::Int)
    Δx, Δy = diff_lag.([x, y], lag=lag) # diff_lag in util.jl
    
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

function ang_btw_vec(v1::Array{<:AbstractFloat,2}, v2::Array{<:AbstractFloat,2})
    timept = min(size(v1, 2), size(v2, 2))
    vec_ang_lst = zeros(timept)
    
    for i = 1 : timept
        vec_ang_lst[i] = vec_ang(v1[:, i], v2[:, i])
    end
    
    vec_ang_lst
end

# take in 2 x num_vec array as iput and return 1 x num_vec array of magnitudes
function magnitude_vec(list_v::Array{<:AbstractFloat,2})
    sqrt.(list_v[1, :] .^ 2 .+ list_v[2, :] .^ 2)
end

# If a to b vector, make into b to a vector, save out the magnitude as well
function reverse_vec(v::Array{<:AbstractFloat,2})
    reversed_v = zeros(2, size(v, 2))
    reversed_v[1, :] = - v[1, :]
    reversed_v[2, :] = - v[2, :]
    
    reversed_v
end

# defining the worm movement vector, by using stage coordinate changes
function mov_vec(x, y, lag::Int)
    timept = size(x, 2)
    mov_v = zeros(2, timept - lag)
    mov_v[1, :] = x[1 : end - lag] 
    mov_v[2, :] = y[1 : end - lag]

    mov_v
end

# A = [0.9 0.1; 0.9 0.1] # transition matrix
# B = [Normal(170,10), Normal(13,5)] # observation dist, mean and variance
function cluster(list::Array{<:AbstractFloat,1}, A, B)
    hmm = HMM(A, B)
    viterbi(hmm, list)
end

function make_vec(x::Array{<:AbstractFloat,1}, y::Array{<:AbstractFloat,1})
    convert(Array{Float64,2}, (hcat(x, y))') # can't do typeof(x) since it's 1D
end

## TODO: tried to make it simpler by doing this, but error saying that the value inside acos is > 1

# # v1 and v2 have to be same dimension array
# function test_ang_btw_vec(v1::Array{<:AbstractFloat,2}, v2::Array{<:AbstractFloat,2})
#     dotprd = sum(dot.(v1, v2), dims=1) # calculate the dot product of two vectors
#     ang_btw_v = acos.(dotprd ./ (magnitude_vec(v1) .* magnitude_vec(v2))) # calculate angle in radians
    
#     ang_btw_v
# end

# # cut the smaller dim vector
# function match_dim(v1::Array{<:AbstractFloat,2}, v2::Array{<:AbstractFloat,2})
#     min_len = min(size(v1,2), size(v2,2))
#     v1 = Array((hcat(v1[1, 1: min_len], v1[2, 1: min_len]))')
#     v2 = Array((hcat(v2[1, 1: min_len], v2[2, 1: min_len]))')
    
#     v1, v2
# end


# 
# """
#     impute_stage(x::Array{<:AbstractFloat,2})

# Imputes missing 2D stage data (NaN or missing) with interpolation

# Arguments
# ---------
# * `pos_stage`: (x,y) location, 2 by T array where T is len(time points)
# """
# function impute_stage(pos_stage::Array{<:AbstractFloat,2})
#     x_imp = impute_list(pos_stage[1, :])
#     y_imp = impute_list(pos_stage[2, :])

#     convert(typeof(pos_stage), (hcat(x_imp,y_imp))')
# end
