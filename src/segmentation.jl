function downsample_unet(img)
    @assert(size(img) == (968, 732))
    img = Flaot32.(img[4:963, 4:723])
    (img[1:2:end-1,1:2:end-1] .+ img[2:2:end,2:2:end] .+
        img[1:2:end-1,2:2:end] .+ img[2:2:end,1:2:end-1]) ./ 4
end

function process_unet_input(img)
    UNet2D.standardize(downsample_unet(img))
end