# function standardize(array::Array{<:AbstractFloat,2})
#     μ = mean(array)
#     σ = std(array)
#     (array .- μ) / σ
# end

# function standardize(array::Array{<:AbstractFloat,3}; dims=(2,3))
#     μ = mean(array, dims=dims)
#     σ = std(array, dims=dims)
#     (array .- μ) ./ σ 
# end

# reshape_array(array::Array{<:AbstractFloat,2}) = reshape(array, (1,1,size(array)...))
# reshape_array(array::Array{<:AbstractFloat,3}) = reshape(array, (size(array,1),1,size(array)[2:3]...))

# function eval_unet(img::Union{Array{Float32,2}, Array{Float32,3}}, model)
#     img_x = reshape_array(img)
#     img_x = torch.from_numpy(img_x).to(device)
    
#     @pywith torch.no_grad() begin
#         img_y_pred = model(img_x)
#         img_y_pred = torch.nn.functional.sigmoid(img_y_pred)
#         img_y_pred = img_y_pred.cpu().numpy()
        
#         if size(img_y_pred)[1] == 1
#             return img_y_pred[1,1,:,:]
#         else
#             return img_y_pred[:,1,:,:]
#         end
#     end
# end

# function create_model(device, path_weights=nothing, n_ch_input=1, n_class=1, n_feature_init=32)
#     # create model
#     model = unet2d.UNet2D(n_ch_input, n_class, n_feature_init, false)
#     !isnothing(path_weights) && model.load_state_dict(torch.load(path_weights,
#             map_location=device))
#     model.eval()
#     model.to(device)
    
#     return model
# end

# function unet(img, threshold)
    
#     if size(img) == (968, 732) # TODO will we be seeing other dimensions of images?
#         img = rotr90(img)
#     end

#     if size(img) == (732, 968)
#         SIZE_X, SIZE_Y = 360, 480; # downsample initial image by half and crop
#         half_x, half_y = 366, 484;
#         diff_x = abs(half_x - SIZE_X)
#         diff_y = abs(half_y - SIZE_Y)

#         resized_img = imresize(img, (366, 484)) # downsize here
#         cropped_img = resized_img[diff_x+1:end, diff_y+1:end]

#         img_float = convert(Array{Float64,2}, cropped_img)

#         img_test_single = Float32.(standardize(img_float))

#         result = eval_unet(img_test_single, model)
#         result_binary = result .> threshold

#         return result, result_binary, img_float
#     end
# end

# function overlay_images(image::Array{Float64,2}, segmentation::BitArray{2}, skeleton::Array{Bool,2}, img_alpha, seg_alpha, skl_alpha)
#     img = np.float32(image)
#     seg = np.float32(segmentation)
#     skl = np.float32(skeleton)
    
#     added_image = cv2.addWeighted(img, img_alpha, seg, seg_alpha, 0)
#     final_image = cv2.addWeighted(added_image, 1, skl, skl_alpha, 0)
    
#     return final_image
# end

# function save_images(image::Array{Float32,2}, directory::String, index::Int64)
#     saving_dir = joinpath(directory, "$(lpad(string(index), 4, "0")).png")
#     plt.imsave(saving_dir, image)
# end

# function encode_movie(input, output; fps=30)
#     run(`ffmpeg -hide_banner -loglevel panic -y -framerate $fps -i $input -c:v libx264 -pix_fmt yuv420p -preset slow -b:v 16M $output`)
#     nothing
# end;

