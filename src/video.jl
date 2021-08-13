"""
    write_behavior_video(path_h5, path_vid=nothing; fps=20, encoder_options=nothing, downsample=true)
Writes a video of the behavior data
Arguments
---------
* `path_h5`: path of the HDF5 file  
* `path_vid`: path of the video to be generated. if `nothing` (default), automatically generates it  
* `fps`: frame rate of the video (default: 20 to match the acquisition rate)  
* `downsample`: if true downsample by factor of 2   
* `encoder_options`: in named tuple (e.g. `(crf="23", preset="slow")`)
### Notes on the encoder options  
`preset`: possible options are: ultrafast, superfast, veryfast, faster, fast, medium – default preset, slow, slower, veryslow.  
`crf`: The range of the CRF scale is 0–51, where 0 is lossless, 23 is the default, and 51 is worst quality possible.  
A lower value generally leads to higher quality, and a subjectively sane range is 17–28.  
Consider 17 or 18 to be visually lossless or nearly so; it should look the same or nearly the same as the input but it isn't technically lossless.  
"""
function write_behavior_video(path_h5, path_vid=nothing; fps=20, encoder_options=nothing, downsample=true)
    path_vid = isnothing(path_vid) ? splitext(path_h5)[1] * "_$(fps)fps.mp4" : path_vid
    if splitext(path_vid)[2] !== ".mp4"
        error("`path_vid` extension must be .mp4")
    end
    
    encoder_options = (crf="23", preset="slow")
    target_pix_fmt = VideoIO.AV_PIX_FMT_YUV420P
    
    h5open(path_h5, "r") do h5f
        img_nir = h5f["img_nir"]
        n_x, n_y, n_t = size(img_nir)

        img_t1 = uint8_to_rgb(downsample ? ds(img_nir[:,:,1]) : img_nir[:,:,1])'
        
        open_video_out(path_vid, img_t1, framerate=20,
            encoder_options=encoder_options, codec_name="libx264",
            target_pix_fmt=target_pix_fmt) do vidf
            @showprogress for t = 2:n_t
                img_ = uint8_to_rgb(downsample ? ds(img_nir[:,:,t]) : img_nir[:,:,t])'
            write(vidf, img_)
            end # t
        end # vid
    end #h5open
end

function uint8_to_rgba(img, alpha=1.0)
    img_ = reinterpret.(N0f8, img)
    RGBA.(img_, img_, img_, fill(N0f8(alpha), size(img_)...))
end

function uint8_to_rgb(img)
    img_ = reinterpret.(N0f8, img)
    RGB.(img_, img_, img_)
end

function ds(img)
    img = Float64.(img)
    round.(UInt8, (img[1:2:end-1,1:2:end-1] .+ img[2:2:end,2:2:end] .+
            img[1:2:end-1,2:2:end] .+ img[2:2:end,1:2:end-1]) ./ 4)
end


function encode_movie(input, output; fps=30)
    run(`ffmpeg -hide_banner -loglevel panic -y -framerate $fps -i $input -c:v libx264 -pix_fmt yuv420p -preset slow -b:v 16M $output`)
    nothing
end;
