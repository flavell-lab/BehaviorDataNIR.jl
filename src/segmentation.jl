euclidean_dist(x1, x2) = norm(x1 .- x2, 2)

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


function get_nn_graph(xs, ys)
    pts_array = hcat(xs, ys)
    nn = py_skl_neighbors.NearestNeighbors(n_neighbors=2).fit(pts_array)
    g = nn.kneighbors_graph()
    g_mat = py_nx.from_scipy_sparse_matrix(g)
    
    g, g_mat
end

function order_pts(xs, ys)
    pts_array = hcat(xs, ys)
    g, g_mat = get_nn_graph(xs, ys)
    list_ord = [collect(py_nx.dfs_preorder_nodes(g_mat, i)) for i = 0:length(xs)-1]
    list_cost = zeros(length(xs))

    for i = 1:length(xs)
        list_cost[i] = sum((pts_array[list_ord[i][1:end-1] .+ 1, :] .- 
                pts_array[list_ord[i][2:end] .+ 1, :]) .^ 2)
    end
    
    list_ord[findmin(list_cost)[2]] .+ 1
end

function longest_shortest(xs, ys)
    g, g_mat = get_nn_graph(xs, ys)
    
    # find terminal nodes
    idx_deg1 = findall(dropdims(sum(g.toarray(), dims=1), dims=1) .!= 2)
    
    # find longest shortest path for end nodes combinations
    list_paths = []
    for nodes = combinations(idx_deg1, 2)
        if py_nx.has_path(g_mat, nodes[1]-1, nodes[2]-1)
            shortest_path = py_nx.shortest_path(g_mat, nodes[1]-1, nodes[2]-1)
            push!(list_paths, (length(shortest_path), shortest_path))
        else
            push!(list_paths, (0, [0]))
        end
    end
    idx_long_short = findmax(map(x->x[1], list_paths))[2]
    
    list_paths[idx_long_short][2] .+ 1 # path_longest_shortest
end

function medial_axis(img_bin, pts_n)
    # medial axis extraction
    img_med_axis = py_ski_morphology.medial_axis(img_bin)
    array_pts = cat(map(x->[x[2], x[1]], findall(img_med_axis))..., dims=2)
    xs = array_pts[2,:]
    ys = array_pts[1,:]

    # find the longest-shortest path
    pts_order = longest_shortest(xs, ys)

    # reorder points
    xs = xs[pts_order]
    ys = ys[pts_order]
    
    # find head/tail and flip
    dist_nose_1 = euclidean_dist((xs[1], ys[1]), pts_n)
    dist_nose_end = euclidean_dist((xs[end], ys[end]), pts_n)

    if dist_nose_1 > dist_nose_end
        reverse!(xs)
        reverse!(ys)
    end

    s = size(img_bin)

    # crop out points that hit the edge of frame
    crop_max = length(xs)
    for i=1:length(xs)
        if xs[i] == s[1] || ys[i] == s[2]
            crop_max = i
            break
        end
    end

    # crop out points in front of the nose
    crop_min = 1
    min_dist = Inf
    for i=1:length(xs)
        d = euclidean_dist((xs[i], ys[i]), pts_n)
        if d < min_dist
            min_dist = d
            crop_min = i
        end
    end


    xs = xs[crop_min:crop_max]
    ys = ys[crop_min:crop_max]

    xs, ys
end

function nn_dist(t, spl)
    spl_pts = spl(t, 1:2)

    sqrt.(dropdims(sum((spl_pts[2:end, :] .- spl_pts[1:end-1, :]) .^ 2, dims=2), dims=2))
end

function cost_dist!(t, spl, dists)
    t = clamp.(t, 0, 1)
    t[1] = 0.
    t[end] = 1.
    dists .= nn_dist(t, spl)
    
    sum((dists .- mean(dists)) .^ 2)
end

function fit_spline(A::Matrix; t=nothing)
    t = isnothing(t) ? LinRange(0,1,size(A,1)) : t
    Interpolations.scale(interpolate(A, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())),
        t, 1:2)
end

function fit_spline(xs, ys, pts_n; n_subsample=15)
    # subsample data points
    spl_data = cat(xs[1:n_subsample:end], ys[1:n_subsample:end], dims=2)
    if (length(xs) - 1) % n_subsample != 0
        spl_data = vcat(spl_data, [xs[end], ys[end]]')
    end

    spl_data[1,:] .= pts_n

    # fit initial spline
    spl_init = fit_spline(spl_data)

    # optimize to uniformly distribute the pts
    t0 = collect(0:0.01:1)
    dists = zeros(Float64, length(t0)-1)
    res_dist = optimize(t -> cost_dist!(t, spl_init, dists), t0, ConjugateGradient(),
        Optim.Options(g_abstol=1e-4))
    x_opt = res_dist.minimizer

    return spl_data, fit_spline(spl_init(x_opt, 1:2))
end
