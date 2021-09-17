"""
Gets the tuning of a neuron's activity pattern to a behavior.
Effectively, this function makes a scatterplot of activity vs behavior,
and for each value of behavior, computes what the typical neuron activity is by
doing a convolution with the neuron activity.

# Arguments:
- `activity`: Neuron activity traces
- `behavior`: Behavior in question
- `conv_fn`: Convolution function to apply to the behavior (eg: `x->pdf(Normal(0,pi/4),x)`)
- `angle` (optional, default `nothing`): If set, the behavior is in terms of an angle. In that case,
    the convolution function will be applied to the difference between the behavioral angle and the `angle` parameter.
"""
function get_tuning(activity, behavior, conv_fn; angle=nothing)
    if !isnothing(angle)
        tuning = median(activity, weights(conv_fn.(recenter_angle.(behavior .- angle))))
        uniform_tuning = median(activity)
    else
        tuning = median(activity, weights(conv_fn.(behavior)))
        uniform_tuning = median(activity)
    end
    return tuning - uniform_tuning
end

"""
Model for Type-1 reversal neurons

# Arguments:
- `a`: Slope of activity rise in response to a reversal
- `b`: Decay rate (0 = no decay)
- `velocity`: Worm velocity
- `init_activity` (default 0): Initial neuron activity
- `v_thresh` (default 0): Velocity threshold
- `vec_to_confocal`: If `velocity` is provided in the NIR time frame, a function that maps NIR to confocal time points.
"""
function reversal_neuron_model(a::Real, b::Real, velocity::Array{<:Real}; init_activity::Real=0, v_thresh::Real=0, vec_to_confocal::Function=identity)
    activity = []
    curr_activity = init_activity
    for t=1:length(velocity)
        curr_activity -= a * min(velocity[t], v_thresh) + b * curr_activity
        push!(activity, curr_activity)
    end
    return vec_to_confocal(activity .- mean(activity))
end


"""
Model for the neuron RIM (aka Type-2 reversal neuron)

# Arguments:
- `a`: Slope of activity rise in response to a reversal
- `b`: Decay rate (0 = no decay)
- `velocity`: Worm velocity
- `init_activity` (default 0): Initial neuron activity
- `v_thresh` (default 0): Velocity threshold
- `vec_to_confocal`: If `velocity` is provided in the NIR time frame, a function that maps NIR to confocal time points.
"""
function RIM_model(a::Real, b::Real, λ::Real, velocity::Array{<:Real}; init_activity::Real=0, v_thresh::Real=0, vec_to_confocal::Function=identity)
    activity = []
    curr_activity = init_activity
    max_t = length(velocity)
    v_avg = ewma(λ, velocity)
            
    for t=1:max_t
        curr_activity -= a * min(v_avg[t], v_thresh) + b * curr_activity
        push!(activity, curr_activity)
    end
    return activity .- mean(activity)
end

"""
Compute exponentially-averaged variable.

# Arguments:
- `λ`: Inverse timescale parameter. 0 = average all previous velocities, Inf = only current timepoint
- `var`: Variable
- `vec_to_confocal`: If `var` is provided in the NIR time frame, a function that maps NIR to confocal time points.
"""
function ewma(λ::Real, var::Array{<:Real}; vec_to_confocal::Function=identity)
    max_t = length(var)
    v_avg = zeros(max_t)
    s = sum([exp(-(max_t-t)*λ) for t=1:max_t])
    v_avg[1] = var[1]
    for t=2:max_t
        v_avg[t] = v_avg[t-1] * (s-1) / s + var[t] / s
    end
    return vec_to_confocal(v_avg)
end

"""
Model for forward neurons

# Arguments:
- `a`: Scaling factor
- `velocity`: Worm's velocity
"""
function forward_neuron_model(a, velocity)
    activity = a .* max.(velocity, 0)
    return activity .- mean(activity)
end

"""
Model for turning and head-angle neurons

# Arguments:
- `a`: Scaling factor
- `b`: Position on the spline to measure angle
- `c`: Reference position on the spline
- `body_angle`: Spline points across all time points
"""
function turning_neuron_model(a, b, c, body_angle)
    return a .* (body_angle[b,:] .- body_angle[c,:])
end
