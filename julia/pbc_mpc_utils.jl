remove(a, b) = a[findall(!in(b), a)]
findin(a, b) = a[findall(in(b), a)]


function get_nodes_by_tag(model::DiscreteModel, tag::String)::Vector{Int64}

    
    labels = get_face_labeling(model)
    tag_id = get_tag_from_name(labels, tag);
    
    dim = 0 # we will look for nodes which has dimension of 0
    
    dface_to_entity = get_face_entity(labels, dim) 
    
    dface_to_isontag = BitVector(undef, num_faces(labels,dim))

    tag_entities = get_tag_entities(labels, tag)

    for i in eachindex(dface_to_entity)
        buf = false
        for entity in tag_entities
            buf += dface_to_entity[i] == entity
        end
        dface_to_isontag[i] = buf
    end

    return findall(dface_to_isontag)

end


# labels = get_face_labeling(model)
# dface_to_entity = get_face_entity(labels, 0)
# tag_entities = get_tag_entities(labels, "back_face")
# dface_to_isontag = BitVector(undef, num_faces(labels,0))
# for entity in tag_entities
#     for i in eachindex(dface_to_entity)
#         dface_to_isontag[i] =  dface_to_entity[i] == entity
#     end
# end
# findall(dface_to_isontag)

function get_node_dofs(fe_space::UnconstrainedFESpace, node_ind)
    triangulation = get_triangulation(fe_space)
    number_of_nodes = num_nodes(triangulation)
    @assert all(node_ind .<= number_of_nodes)
    
    dims = num_point_dims(triangulation)

    # ndofs = num_free_dofs(fe_space);

    dofs = Vector{Int64}(undef, length(node_ind)*dims);

    for i in eachindex(node_ind)
        for j in 1:dims
            dofs[i*dims - 3 + j] = node_ind[i] * dims - 3 + j
        end
    end
    return dofs
end


"""
function applyControlDependeMPC!(  
    mK::AbstractMatrix{<:Number}, 
    vF::AbstractVector{<:Number}, 
    primary_dofs::AbstractVector{<:Integer}, 
    secondary_dofs::AbstractVector{<:Integer}, 
    koefs::AbstractVector{<:Number}, 
    offset::AbstractVector{<:Number}
)
    # sec[i] = prim[i] * koef[i] + offs[i]
    # 1 * s - k * p = o
    # 
    # | 1  -k | =  1*o
    # |-k  k^2| = -k*o
    @assert size(primary_dofs, 1) == size(secondary_dofs,1)
    @assert size(primary_dofs, 1) == size(koefs,  1)
    @assert size(primary_dofs, 1) == size(offset, 1)


    for el in eachindex(primary_dofs)
        i = secondary_dofs[el];
        j = primary_dofs[el];
        mW = [
            1 -koeffs[el];
            -koeffs[el] koeffs[el]^2
        ] * ω;

        mK[[i,j],[i,j]] += mW;
        vF[[i,j]] += offset[el] * [1 , -koeffs[el]] * ω;
    end

end




function add_relation!(mpc_map, master, slave, coeff, offset)
    if haskey(mpc_map, master)
        uber_master, uber_coeff, uber_offset = mpc_map[master]
        mpc_map[slave] = (uber_master, coeff * uber_coeff, coeff * uber_offset + offset)
    elseif haskey(mpc_map, slave)
        println("attention!\n   was: \$(slave) -> \$(mpc_map[slave])")
        mpc_map[slave] = (master, coeff, offset)
        println("became: \$(slave) -> \$(mpc_map[slave])")
    else
        mpc_map[slave] = (master, coeff, offset)
    end
end
"""

function build_relation_map(tag1::String, tag2::String, model::DiscreteModel, relation::Function)
    first_nodes = get_nodes_by_tag(model, tag1);
    second_nodes = get_nodes_by_tag(model, tag2);
    
    @assert length(first_nodes) == length(second_nodes) "Size of node sets doesn't match"

    node_coordinates = get_node_coordinates(model);

    first_nodes_sorted = sortperm(node_coordinates[first_nodes]);
    second_nodes_sorted = sortperm(node_coordinates[second_nodes]);

    check = BitVector(undef, length(first_nodes))
    # ch = Vector{Float64}(undef, length(first_nodes))
    for (i, (first_ind, second_ind)) in enumerate(zip(first_nodes[first_nodes_sorted], second_nodes[second_nodes_sorted]))
        check[i] = relation(node_coordinates[first_ind], node_coordinates[second_ind]) < 1e-3
    end

    # return ch
    
    @assert all(check) "Node sets don't match the relationship"

    return (first_nodes[first_nodes_sorted], second_nodes[second_nodes_sorted])
end

function build_relation_map(tags1::Tuple{String, String, String}, tags2::Tuple{String, String, String}, model::DiscreteModel, relation::Function)
    first_nodes = Vector{Int64}()
    second_nodes = Vector{Int64}()

    for tag1 in tags1
        push!(first_nodes, get_nodes_by_tag(model, tag1)...);
    end

    for tag2 in tags2
        push!(second_nodes, get_nodes_by_tag(model, tag2)...);
    end


    @assert length(first_nodes) == length(second_nodes) "Size of node sets doesn't match"

    node_coordinates = get_node_coordinates(model);

    first_nodes_sorted = sortperm(node_coordinates[first_nodes]);
    second_nodes_sorted = sortperm(node_coordinates[second_nodes]);

    check = BitVector(undef, length(first_nodes))
    # ch = Vector{Float64}(undef, length(first_nodes))
    for (i, (first_ind, second_ind)) in enumerate(zip(first_nodes[first_nodes_sorted], second_nodes[second_nodes_sorted]))
        check[i] = relation(node_coordinates[first_ind], node_coordinates[second_ind]) < 1e-5
    end

    # return ch
    
    @assert all(check) "Node sets don't match the relationship"

    return (first_nodes[first_nodes_sorted], second_nodes[second_nodes_sorted])
end