function create_plate_selection_bbox(dir, pos, size, eps = 1e-3)
    
    mask = [1,2,3] .== dir;
    bbox_start = Vector{Float64}(undef, 3);
    bbox_end = Vector{Float64}(undef, 3);

    bbox_start[mask] .= pos;
    bbox_start[.!mask] .= -size/2;
    bbox_start .-= eps;
    
    bbox_end[mask] .= pos;
    bbox_end[.!mask] .= size/2;
    bbox_end .+= eps;

    bbox = [
        bbox_start;
        bbox_end;
    ]
    return bbox;
end



function prepare_to_periodic_sphere(xc, yc, zc, r, cs)
    centers = [[xc, yc, zc]];
       
    if abs(cs/2 - xc) < r
        push!(centers, [xc - cs, yc, zc])
    elseif abs(-cs/2 - xc) < r
        push!(centers, [xc + cs, yc, zc])
    end
    n = length(centers)
    if abs(cs/2 - yc) < r
        for i in 1:n
            pc = centers[i]
            push!(centers, [pc[1], pc[2] - cs, pc[3]])
        end
    elseif abs(-cs/2 - yc) < r
        for i in 1:n
            pc = centers[i]
            push!(centers, [pc[1], pc[2] + cs, pc[3]])
        end
    end
    n = length(centers)
    if abs(cs/2 - zc) < r
        for i in 1:n
            pc = centers[i]
            push!(centers, [pc[1], pc[2],pc[3] - cs])
        end
    elseif abs(-cs/2 - zc) < r
        for i in 1:n
            pc = centers[i]
            push!(centers, [pc[1], pc[2], pc[3] + cs])
        end
    end
    return centers
end

function prepare_to_periodic(xc, yc, zc, r, cs)
    tr_list = zeros(3)
    if abs(cs/2 - xc) < r
        tr_list[1] = -1
    elseif abs(-cs/2 - xc) < r
        tr_list[1] = 1
    end
    if abs(cs/2 - yc) < r
        tr_list[2] = -1
    elseif abs(-cs/2 - yc) < r
        tr_list[2] = 1
    end
    if abs(cs/2 - zc) < r
        tr_list[3] = -1
    elseif abs(-cs/2 - zc) < r
        tr_list[3] = 1
    end
    return tr_list
end

using Test

@test begin
    pc = [0.0, 0.0, 0.0]
    cs = 1
    r = 0.1
    prepare_to_periodic(pc..., r, cs) == [0,0,0]
end

@test begin
    pc = [0.5, 0.0, 0.0]
    cs = 1
    r = 0.1
    prepare_to_periodic(pc..., r, cs) == [-1,0,0]
end

@test begin
    pc = [0.5, -0.5, 0.0]
    cs = 1
    r = 0.1
    prepare_to_periodic(pc..., r, cs) == [-1,1,0]
end

function compute_transformation_matrices_for_periodic(list, cs)
    # [1 0 0] - Tx
    # [0 1 0] - Ty
    # [0 0 1] - Tz

    # [1 1 0] - Tx*Ty
    # [0 1 1] - Ty*Tz
    # [1 0 1] - Tx*Tz

    # [1 1 1] - Tx*Ty*Tz

    tr_num = Int64(2^sum(abs.(list)) - 1)
    tr_matrices = Array{Float64, 3}(undef, 4, 4, tr_num)
    I = Float64[
        1 0 0 0;
        0 1 0 0;
        0 0 1 0;
        0 0 0 1
    ]
    k = 0;
    Tx = copy(I); Tx[1, 4] = list[1] * cs;
    Ty = copy(I); Ty[2, 4] = list[2] * cs;
    Tz = copy(I); Tz[3, 4] = list[3] * cs;
    if list[1] != 0
        k+=1
        tr_matrices[:,:,k] = Tx
    end
    if list[2] != 0
        k+=1
        tr_matrices[:,:,k] = Ty
    end
    if list[3] != 0
        k+=1
        tr_matrices[:,:,k] = Tz
    end
    if all(list[[1,2]] .!= 0)
        k+=1
        tr_matrices[:,:,k] = Tx*Ty
    end
    if all(list[[2,3]] .!= 0)
        k+=1
        tr_matrices[:,:,k] = Ty*Tz
    end
    if all(list[[1,3]] .!= 0)
        k+=1
        tr_matrices[:,:,k] = Tx*Tz
    end
    if all(list .!= 0)
        k+=1
        tr_matrices[:,:,k] = Tx*Ty*Tz
    end
    return tr_matrices

end
begin
    pc = [0.5, -0.5, 0.5]
    cs = 1
    r = 0.1
    list = prepare_to_periodic(pc..., r, cs)
    matrices = prepare_transformation_list_for_periodic(list, 1)

    cmp = matrices[:,:,7]
    reshape(cmp', 16)'
    
end




function add_periodical_sphere(xc, yc, zc, r, cs)
    
    list = prepare_to_periodic(xc, yc, zc, r, cs)
    tr_matrices = compute_transformation_matrices_for_periodic(list, cs)
    n = size(tr_matrices, 3);
    
    pores = Vector{Tuple{Int64, Int64}}()
    root_pore = model.occ.add_sphere(xc, yc, zc, r)
    push!(pores, (3, root_pore))
    for i in 1:n
        pore_id = model.occ.copy((3, root_pore))[1][1]
        model.occ.affine_transform((3, root_pore), reshape(tr_matrices[:,:, i]',16))
        push!(pores, (3,pore_id))
    end
    return pores
end

function make_body_periodical(tag, cs)
    body_bb = model.occ.getBoundingBox(3, tag)
    
end

