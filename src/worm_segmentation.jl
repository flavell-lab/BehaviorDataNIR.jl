function standardize(array::Array{<:AbstractFloat,2})
    μ = mean(array)
    σ = std(array)
    (array .- μ) / σ
end

function standardize(array::Array{<:AbstractFloat,3}; dims=(2,3))
    μ = mean(array, dims=dims)
    σ = std(array, dims=dims)
    (array .- μ) ./ σ 
end
