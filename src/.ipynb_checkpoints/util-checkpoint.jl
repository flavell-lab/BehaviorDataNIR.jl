function diff_lag(A::AbstractVector; lag::Int)
    A[lag+1:end] .- A[1:end-lag]
end