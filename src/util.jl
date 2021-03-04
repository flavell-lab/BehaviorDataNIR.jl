function diff_lag(A::AbstractVector, lag::Int)
    A[lag+1:end] .- A[1:end-lag]
end

function vec_ang(v1, v2)
    acos(dot(v1, v2) / (norm(v1, 2) * norm(v2, 2)))
end