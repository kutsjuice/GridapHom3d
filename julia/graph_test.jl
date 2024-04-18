using Graphs

left_nodes = [1,2,3]
right_nodes = [7,6,5]

top_nodes = [1,8,7]
bottom_nodes = [3,4,5]


left_dofs = [1,2,3,4,5,6]
right_dofs = [13,14, 11, 12, 9, 10]

top_dofs = [1, 2, 15, 16, 13, 14]
bottom_dofs = [5, 6, 7, 8, 9, 10]

findin(a, b) = a[findall(in(b), a)]


lt_dofs = findin(left_dofs, top_dofs)
lb_dofs = findin(left_dofs, bottom_dofs)

cell_size = 1;
domain = (-cell_size/2, cell_size/2, -cell_size/2, cell_size/2, -cell_size/2, cell_size/2)
partition = (2, 2, 2)
using Gridap

model = CartesianDiscreteModel(domain,partition)
vtk_file = "created_mesh" |> input;
writevtk(model, vtk_file);


back_tag = "tag_21"
front_tag = "tag_22"
bottom_tag = "tag_23"
top_tag = "tag_24"
left_tag = "tag_25"
right_tag = "tag_26"


include("pbc_mpc_utils.jl")

back_selector(point::VectorValue) = abs(point[3] + cell_size/2) < 1e-6

front_selector(point::VectorValue) = abs(point[3] - cell_size/2) < 1e-6

left_selector(point::VectorValue) = abs(point[1] + cell_size/2) < 1e-6

right_selector(point::VectorValue) = abs(point[1] - cell_size/2) < 1e-6

bottom_selector(point::VectorValue) = abs(point[2] + cell_size/2) < 1e-6

top_selector(point::VectorValue) = abs(point[2] - cell_size/2) < 1e-6

function select_nodes_by_selector(model::DiscreteModel, selector)
    return findall(lazy_map(selector, get_node_coordinates(model)))    
end


function select_nodes_by_selector(coords::Vector{VectorValue{3, Float64}}, selector::Function)
    return findall(lazy_map(selector, coords))    
end

coords = collect(reshape(get_node_coordinates(model), 27))

calc_node_dofs_3D(ind::Integer) = (ind*3 -2, ind*3 -1, ind*3)

function calc_node_dofs_3D(nodes::Vector{<:Integer})
    dofs = Vector{Int64}(undef, 3*length(nodes))
    for i in eachindex(nodes)
        dofs[3i - 2] = nodes[i] * 3 - 2;
        dofs[3i - 1] = nodes[i] * 3 - 1;
        dofs[3i - 0] = nodes[i] * 3 - 0;
    end
    return dofs
end


lft_face_nodes = select_nodes_by_selector(coords, left_selector)
rht_face_nodes = select_nodes_by_selector(coords, right_selector)
btm_face_nodes = select_nodes_by_selector(coords, bottom_selector)
top_face_nodes = select_nodes_by_selector(coords, top_selector)
bck_face_nodes = select_nodes_by_selector(coords, back_selector)
frt_face_nodes = select_nodes_by_selector(coords, front_selector)



# In overal, we hav 12 edges:
#  1. left-top
#  2. left-bottom
#  3. left-front
#  4. left-back
#  5. right-top
#  6. right-bottom
#  7. right-front
#  8. right-back
#  9. front-top
# 10. front-bottom
# 11. back-top
# 12. back-bottom
#
# and 8 vertices
# 1. left-top-front
# 2. left-top-back
# 3. left-bottom-front
# 4. left-bottom-back
# 5. right-top-front
# 6. right-top-back
# 7. right-bottom-front
# 8. right-bottom-back


lft_top_edge_nodes = findin(lft_nodes, top_nodes)
lft_btm_edge_nodes = findin(lft_nodes, btm_nodes)
lft_frt_edge_nodes = findin(lft_nodes, frt_nodes)
lft_bck_edge_nodes = findin(lft_nodes, bck_nodes)
rht_top_edge_nodes = findin(rht_nodes, top_nodes)
rht_btm_edge_nodes = findin(rht_nodes, btm_nodes)
rht_frt_edge_nodes = findin(rht_nodes, frt_nodes)
rht_bck_edge_nodes = findin(rht_nodes, bck_nodes)
top_frt_edge_nodes = findin(top_nodes, frt_nodes)
top_bck_edge_nodes = findin(top_nodes, bck_nodes)
btm_frt_edge_nodes = findin(btm_nodes, frt_nodes)
btm_bck_edge_nodes = findin(btm_nodes, bck_nodes)


lft_top_frt_vertex_node = findin(lft_top_edge_nodes, frt_face_nodes)
lft_top_bck_vertex_node = findin(lft_top_edge_nodes, bck_face_nodes)
lft_btm_frt_vertex_node = findin(lft_btm_edge_nodes, frt_face_nodes)
lft_btm_bck_vertex_node = findin(lft_btm_edge_nodes, bck_face_nodes)

rht_top_frt_vertex_node = findin(rht_top_edge_nodes, frt_face_nodes)
rht_top_bck_vertex_node = findin(rht_top_edge_nodes, bck_face_nodes)
rht_btm_frt_vertex_node = findin(rht_btm_edge_nodes, frt_face_nodes)
rht_btm_bck_vertex_node = findin(rht_btm_edge_nodes, bck_face_nodes)

# next - remove edges from faces and then - vertices from edges

lft_face_nodes = remove(lft_face_nodes, lft_bck_edge_nodes)
lft_face_nodes = remove(lft_face_nodes, lft_frt_edge_nodes)
lft_face_nodes = remove(lft_face_nodes, lft_top_edge_nodes)
lft_face_nodes = remove(lft_face_nodes, lft_btm_edge_nodes)

rht_face_nodes = remove(rht_face_nodes, rht_bck_edge_nodes)
rht_face_nodes = remove(rht_face_nodes, rht_frt_edge_nodes)
rht_face_nodes = remove(rht_face_nodes, rht_top_edge_nodes)
rht_face_nodes = remove(rht_face_nodes, rht_btm_edge_nodes)

btm_face_nodes = remove(btm_face_nodes, lft_btm_edge_nodes)
btm_face_nodes = remove(btm_face_nodes, rht_btm_edge_nodes)
btm_face_nodes = remove(btm_face_nodes, btm_bck_edge_nodes)
btm_face_nodes = remove(btm_face_nodes, btm_frt_edge_nodes)

top_face_nodes = remove(top_face_nodes, lft_top_edge_nodes)
top_face_nodes = remove(top_face_nodes, rht_top_edge_nodes)
top_face_nodes = remove(top_face_nodes, top_bck_edge_nodes)
top_face_nodes = remove(top_face_nodes, top_frt_edge_nodes)

bck_face_nodes = remove(bck_face_nodes, lft_bck_edge_nodes)
bck_face_nodes = remove(bck_face_nodes, rht_bck_edge_nodes)
bck_face_nodes = remove(bck_face_nodes, btm_bck_edge_nodes)
bck_face_nodes = remove(bck_face_nodes, top_bck_edge_nodes)

frt_face_nodes = remove(frt_face_nodes, lft_frt_edge_nodes)
frt_face_nodes = remove(frt_face_nodes, rht_frt_edge_nodes)
frt_face_nodes = remove(frt_face_nodes, btm_frt_edge_nodes)
frt_face_nodes = remove(frt_face_nodes, top_frt_edge_nodes)


lft_top_edge_nodes = remove(lft_top_edge_nodes, lft_top_bck_vertex_node)
lft_top_edge_nodes = remove(lft_top_edge_nodes, lft_top_frt_vertex_node)
lft_btm_edge_nodes = remove(lft_btm_edge_nodes, lft_btm_bck_vertex_node)
lft_btm_edge_nodes = remove(lft_btm_edge_nodes, lft_btm_frt_vertex_node)
lft_frt_edge_nodes = remove(lft_frt_edge_nodes, lft_top_frt_vertex_node)
lft_frt_edge_nodes = remove(lft_frt_edge_nodes, lft_btm_frt_vertex_node)
lft_bck_edge_nodes = remove(lft_bck_edge_nodes, lft_top_bck_vertex_node)
lft_bck_edge_nodes = remove(lft_bck_edge_nodes, lft_btm_bck_vertex_node)
rht_top_edge_nodes = remove(rht_top_edge_nodes, rht_top_bck_vertex_node)
rht_top_edge_nodes = remove(rht_top_edge_nodes, rht_top_frt_vertex_node)
rht_btm_edge_nodes = remove(rht_btm_edge_nodes, rht_btm_bck_vertex_node)
rht_btm_edge_nodes = remove(rht_btm_edge_nodes, rht_btm_frt_vertex_node)
rht_frt_edge_nodes = remove(rht_frt_edge_nodes, rht_top_frt_vertex_node)
rht_frt_edge_nodes = remove(rht_frt_edge_nodes, rht_btm_frt_vertex_node)
rht_bck_edge_nodes = remove(rht_bck_edge_nodes, rht_top_bck_vertex_node)
rht_bck_edge_nodes = remove(rht_bck_edge_nodes, rht_btm_bck_vertex_node)
top_frt_edge_nodes = remove(top_frt_edge_nodes, lft_top_frt_vertex_node)
top_frt_edge_nodes = remove(top_frt_edge_nodes, rht_top_frt_vertex_node)
top_bck_edge_nodes = remove(top_bck_edge_nodes, lft_top_bck_vertex_node)
top_bck_edge_nodes = remove(top_bck_edge_nodes, rht_top_bck_vertex_node)
btm_frt_edge_nodes = remove(btm_frt_edge_nodes, lft_btm_frt_vertex_node)
btm_frt_edge_nodes = remove(btm_frt_edge_nodes, rht_btm_frt_vertex_node)
btm_bck_edge_nodes = remove(btm_bck_edge_nodes, lft_btm_bck_vertex_node)
btm_bck_edge_nodes = remove(btm_bck_edge_nodes, rht_btm_bck_vertex_node)


# create maps; in overal we should have 3 facet maps + 9 edge maps + 7  vertex maps
# and for each  we should have relation func

btm_face_2_top_face_relation(btm_node::VectorValue, top_node::VectorValue) = norm(top_node - btm_node - VectorValue(0,0,cell_size)) 




pbc_mpc = Dict{Int64, Tuple{Int64, Float64, Float64}}()



cords[lft_bck_edge_nodes]



# 1. add all relations for cell vertices nodes - all nodes are dependent from uber-uber-master node
# 2. add relations for edge nodes - select 3 sets of uber

add_relation!(pbc_mpc, )