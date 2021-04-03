function downsample_unet_input(img)
    @assert(size(img) == (968, 732))
    img = Float32.(img[4:963, 4:723])
    (img[1:2:end-1,1:2:end-1] .+ img[2:2:end,2:2:end] .+
        img[1:2:end-1,2:2:end] .+ img[2:2:end,1:2:end-1]) ./ 4
end

function largest_segment!(img_bin, img_label)
    img_label .= 0
    label_2d!(img_bin, img_label)
    list_seg = get_segmented_instance(img_label, img_bin, img_bin)
    largest_obj_id = sort(list_seg, by=x->x.area, rev=true)[1].obj_id
    
    img_label .== largest_obj_id
end

"""
    segment_worm!(model, img_raw, img_label; θ=0.75)

Segment the worm given img (size: (968, 732)).

Arguments
---------
* `model`: unet model
* `img_raw`: raw img
* `img_label`: pre-allocated label array <:Int32. Size: (480,360)
Optional arguments
--------
* `θ=0.75`: unet output threshold
"""
function segment_worm!(model, img_raw, img_label; θ=0.75)
    @assert(eltype(img_label) == Int32)
    @assert(size(img_label) == (480,360))
    img_raw_ds = downsample_unet_input(img_raw)
    img_unet = eval_model(UNet2D.standardize(img_raw_ds), model)
    img_bin = closing(img_unet .> θ)
    
    img_raw_ds, largest_segment!(img_bin, img_label)
end