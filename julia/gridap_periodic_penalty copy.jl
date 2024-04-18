using Gridap
using GridapGmsh
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.ReferenceFEs
# using Gmsh
using Gridap.Visualization

# using WriteVTK

using LinearAlgebra
# using SparseArrays
# using LinearSolve

include("pbc_mpc_utils.jl")
include("penalty_mpc.jl")

input(filename::String) = "geometry/" * filename;

msh_file = "periodic.msh" |> input
# gmsh.initialize()
model = GmshDiscreteModel(msh_file, renumber=false);

vtk_file = "imported_mesh" |> input;
writevtk(model, vtk_file);
##


# get labels

# labels = get_face_labeling(model)

# dim = 0 # we will look for nodes which has dimension of 0
# dface_to_entity = get_face_entity(labels, dim) # wtf

# # cdat = []
# tag = 1 # i suppose, we will lok for "center" tags
# get_tag_name(labels, tag)

# dface_to_isontag = BitVector(undef, num_faces(labels,dim))

# tag_entities = get_tag_entities(labels, tag)

# for entity in tag_entities
#     for i in eachindex(dface_to_entity)
#         dface_to_isontag[i] =  dface_to_entity[i] == entity
#     end
# end
# findall(dface_to_isontag)



# get_nodes_by_tag(model, "cell")

##

const E = 70.0e9
const ν = 0.33
const λ = (E*ν)/((1+ν)*(1-2*ν))
const μ = E/(2*(1+ν))


σ(ε) = λ*tr(ε)*one(ε) + 2*μ*ε

order = 1
degree = order*2
Ω = Triangulation(model)


dΩ = Measure(Ω,degree)

reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},1)
test_space_V = TestFESpace(model,reffe,conformity=:H1)#,dirichlet_tags = ["center"])

get_nodes_by_tag(model, "left_face")

##
dirichelt_nodes = get_nodes_by_tag(model, "dirichelet_tag")
a(u, v) = 
l(v) = 
op = AffineFEOperator(a, l, test_space, trial_space)
mK = get_matrix(op)

dirichelet_dofs = calc_dofs_for_nodes_3D(dirichelt_nodes)

1 -> 1, 2, 3

i -> 3i-2, 3i-1, 3i


get_cell_dof_ids -> Vector из 12 компонент
[
    x_dof_1
    x_dof_2
    x_dof_3
    x_dof_4
    y_dof_1
    y_dof_2
    y_dof_3
    y_dof_4
    z_dof_1
    z_dof_2
    z_dof_3
    z_dof_4
]

get_cell_node_ids -> Vector 4
[
    node_1
    node_2
    node_3
    node_4
]


body1 -> mK1
body2 -> mK2

mK = [
    mK1 0;
    0   mK2;
]

узлы тела 1 - 1:N
узлы тела 2 - 1:M

узлы сборки - 1:(N+M)
узлы тела 1 - 1:N
узлы тела 2 - N+1:N+M


coordinates_body1 = get_node_coordnates(model1)
coordinates_body2 = get_node_coordnates(model2)


connection_face_body1_nodes = [13, 45, 56]
connection_face_body2_nodes = [13, 14, 15] .+ N

coordinates_body1[connection_face_body1_nodes] - координаты узлов в соединении
coordinates_body2[connection_face_body2_nodes] - координаты узлов в соединении

connection_map = [
    13 -> 14
    45 -> 13
    56 -> 15
]




6 CЧ = 0;
6 CФ - формы rigid body

##

# cell_nodes = get_cell_node_ids(Ω)
# cell_dofs = get_cell_dof_ids(test_space_V)

# cell_dofs[1]
# cell_nodes[1]




trial_space_U = TrialFESpace(test_space_V)

##
# weak form 

bilinear_form_stiffness(u,v) = ∫( ε(v) ⊙ (σ∘ε(u)) )*dΩ
linear_form_stiffness(v) = 0; #∫(dot(f, v))*dΩ


operator_K = AffineFEOperator(bilinear_form_stiffness, linear_form_stiffness, trial_space_U, test_space_V)




stiffness_matrix = get_matrix(operator_K)
rhs_vector = get_vector(operator_K)


# writevtk(Ω, res_file ,  cellfields=["uh"=>uh_lin])






function left_to_right_relation(left::VectorValue, right::VectorValue)
    return norm( right - left - VectorValue(0.5, 0.0, 0.0))
end

function back_to_front_relation(back::VectorValue, front::VectorValue)
    return norm(front - back - VectorValue(0.0, 0.0, 0.5))
end

function bottom_to_top_relation(bottom::VectorValue, top::VectorValue)
    return norm(top - bottom - VectorValue(0.0, 0.5, 0.0))
end

left_to_right_map = build_relation_map("left_face", "right_face", model, left_to_right_relation)
buttom_to_top_map = build_relation_map("bottom_face", "top_face", model, bottom_to_top_relation)
back_to_front_map = build_relation_map("back_face", "front_face", model, back_to_front_relation)

# maximum(left_to_right_map)
# maximum(buttom_to_top_map)
# maximum(back_to_front_map)


function build_xx_tension_pbc_mpc_map(l2r, d2u, b2f, displacement)

end

function build_yz_shear_pbc_mpc_map()

end

#
#  1--------8--------7
#  |       /|\       |
#  |        |        |
#  2 <-------------> 6
#  |        |        |
#  |       \|/       |
#  3--------4--------5
#   u3 = u1 -> 5 = 1
#   u4 = u8 -> 7 = 15
#   u5 = u7 -> 9 = 13
#   u2 = u6 -> 3 = 11
#   u1 = u7 -> 1 = 13
#   -----------------


b2t = [
    3 1;
    4 8;
    5 7;
]

l2r = [
    1 7;
    2 6;
    3 5;
]

mpc_map = Dict{Int, Tuple{Int, Float64, Float64}}()



for row in eachrow(b2t)

    vi = row[1]*2;
    ui = vi - 1;

    vj = row[2]*2;
    uj = vj - 1;

    add_relation!(mpc_map, uj, ui, 1.0, 0.0);
    add_relation!(mpc_map, vj, vi, 1.0, 0.1);
end

for row in eachrow(l2r)
    vi = row[1]*2;
    ui = vi - 1;

    vj = row[2]*2;
    uj = vj - 1;

    add_relation!(mpc_map, uj, ui, 1.0, 0.0);
    add_relation!(mpc_map, vj, vi, 1.0, 0.0);
end
mpc_map







function resolve_dependencies()
end