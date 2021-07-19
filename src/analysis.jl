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
Model for type-1 reversal neurons

# Arguments:
- `a`: Slope of activity rise in response to a reversal
- `b`: Decay rate (0 = instantaneous, 1 = no decay)
- `reversals`: Time points where the worm reverses
- `max_t`: Maximum time point
- `init_activity` (default 0): Initial neuron activity
"""
function reversal_neuron_model_1(a, b, reversals, max_t; init_activity=0)
    activity = []
    curr_activity = init_activity
    for t=1:max_t
        if t in reversals
            curr_activity += a
        else
            curr_activity *= b
        end
        push!(activity, curr_activity)
    end
    return activity .- mean(activity)
end    

"""
Model for type-2 reversal neurons

# Arguments:
- `a`: Slope of activity rise in response to a reversal
- `b`: Decay rate (0 = instantaneous, 1 = no decay)
- `c`: Steepness of neuron activation dependence on reversal duration
- `threshold`: Reversal duration threshold for 1/2 activity
- `reversals`: Time points where the worm reverses
- `max_t`: Maximum time point
- `init_activity` (default 0): Initial neuron activity
"""
function reversal_neuron_model_2(a, b, c, threshold, reversals, max_t; init_activity=0)
    activity = []
    curr_activity = init_activity
    rev_num = 0
    curr_a = a
    for t=1:max_t
        if t in reversals
            if rev_num == 0
                rev_num = minimum([T for T in t:max_t if !(T in reversals)]) - t + 1
                curr_a = a * (c^rev_num) / (c^rev_num + c^threshold)
            end
            curr_activity += curr_a
        else
            curr_activity *= b
            rev_num = 0
        end
        push!(activity, curr_activity)
    end
    return activity .- mean(activity)
end   
