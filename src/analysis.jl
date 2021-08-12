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
Model for reversal neurons

# Arguments:
- `a`: Slope of activity rise in response to a reversal
- `b`: Decay rate (0 = instantaneous, 1 = no decay)
- `len_thresh`: Ignore reversal events with duration at most this. For type-1 neurons, set to 0.
- `rev_lengths`: For each time point, 0 if the worm is going forward, or the duration of that reversal event if it's reversing
- `max_t`: Maximum time point
- `init_activity` (default 0): Initial neuron activity
"""
function reversal_neuron_model(a, b, len_thresh, rev_lengths, max_t; init_activity=0)
    activity = []
    curr_activity = init_activity
    curr_a = a
    for t=1:max_t
        if rev_lengths[t] > len_thresh
            curr_activity += curr_a
        else
            curr_activity *= b
        end
        push!(activity, curr_activity)
    end
    return activity .- mean(activity)
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
