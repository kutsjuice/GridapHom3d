using Gridap
using GridapGmsh
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.Visualization

using LinearAlgebra
using LinearSolve

include("pbc_mpc_utils.jl")
include("penalty_mpc.jl")

input(filename::String) = "geometry/" * filename;
results(filename::String) = "results/" * filename;

msh_file = "periodic.msh" |> input
# gmsh.initialize()
model = GmshDiscreteModel(msh_file, renumber=false)

vtk_file = "imported_mesh" |> input;
writevtk(model, vtk_file);
##
first_nodes = get_nodes_by_tag(model, "left_face")
second_nodes = get_nodes_by_tag(model, "right_face")

@assert length(first_nodes) == length(second_nodes) "Size of node sets doesn't match"

node_coordinates = get_node_coordinates(model)

first_nodes_crd = node_coordinates[first_nodes]
second_nodes_crd = node_coordinates[second_nodes]

first_nodes_sorted = sortperm(node_coordinates[first_nodes])
second_nodes_sorted = sortperm(node_coordinates[second_nodes])
node_coordinates[first_nodes]
node_coordinates[second_nodes]

function get_nodes_in_bounding_box(model, bounding_box)
    node_coordinates = get_node_coordinates(model)
    nodes_in_bbox = BitVector(undef, length(node_coordinates))

    for i in eachindex(nodes_in_bbox)
        nodes_in_bbox[i] = 
            node_coordinates[i][1] >= bounding_box[1] && 
            node_coordinates[i][1] <= bounding_box[4] &&
            
            node_coordinates[i][2] >= bounding_box[2] && 
            node_coordinates[i][2] <= bounding_box[5] &&
            
            node_coordinates[i][3] >= bounding_box[3] && 
            node_coordinates[i][3] <= bounding_box[6];
    end
    return findall(nodes_in_bbox)
end

cell_size = 0.5
left_face_bbox = [
    -cell_size/2 - 1e-3,
    -cell_size/2 - 1e-3,
    -cell_size/2 - 1e-3,

    -cell_size/2 + 1e-3, # <<==
    cell_size/2 + 1e-3,
    cell_size/2 + 1e-3,
]
left_nodes = get_nodes_in_bounding_box(model, left_face_bbox)
coords = get_node_coordinates(model)
coords[left_nodes]
b = [p[1] == -0.25 for p in coords[left_nodes]]
all(b)

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


trial_space_U = TrialFESpace(test_space_V)
##
# weak form 

bilinear_form_stiffness(u,v) = ∫( ε(v) ⊙ (σ∘ε(u)) )*dΩ
linear_form_stiffness(v) = 0; #∫(dot(f, v))*dΩ


operator_K = AffineFEOperator(bilinear_form_stiffness, linear_form_stiffness, trial_space_U, test_space_V)




stiffness_matrix = get_matrix(operator_K)
rhs_vector = get_vector(operator_K)


# writevtk(Ω, res_file ,  cellfields=["uh"=>uh_lin])

##

# rank(stiffness_matrix)

##
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
bottom_to_top_map = build_relation_map("bottom_face", "top_face", model, bottom_to_top_relation)
back_to_front_map = build_relation_map("back_face", "front_face", model, back_to_front_relation)

n1 = length(left_to_right_map[1])
n2 = length(bottom_to_top_map[1])
n3 = length(back_to_front_map[1])

# get_nodes_by_tag(model, "left_face")

i = 1;
dof_map = Dof2DofMap((n1+n2+n3)*3)

(n1+n2+n3)*3
for k in eachindex(left_to_right_map[1])
    # X axis
    dof_map.map[i] = (left_to_right_map[1][k]*3 - 2, left_to_right_map[2][k]*3 - 2);
    dof_map.offsets[i] = 0.1;
    i += 1;
    # Y axis
    dof_map.map[i] = (left_to_right_map[1][k]*3 - 1, left_to_right_map[2][k]*3 - 1);
    dof_map.offsets[i] = 0.0;
    i += 1;
    # Z axis
    dof_map.map[i] = (left_to_right_map[1][k]*3, left_to_right_map[2][k]*3);
    dof_map.offsets[i] = 0.0;
    i += 1;
end
for k in eachindex(bottom_to_top_map[1])
    # X axis
    dof_map.map[i] = (bottom_to_top_map[1][k]*3 - 2, bottom_to_top_map[2][k]*3 - 2);
    dof_map.offsets[i] = 0.0;
    i += 1;
    # Y axis
    dof_map.map[i] = (bottom_to_top_map[1][k]*3 - 1, bottom_to_top_map[2][k]*3 - 1);
    dof_map.offsets[i] = 0.0;
    i += 1;
    # Z axis
    dof_map.map[i] = (bottom_to_top_map[1][k]*3, bottom_to_top_map[2][k]*3);
    dof_map.offsets[i] = 0.0;
    i += 1;
end
for k in eachindex(back_to_front_map[1])
    # X axis
    dof_map.map[i] = (back_to_front_map[1][k]*3 - 2, back_to_front_map[2][k]*3 - 2);
    dof_map.offsets[i] = 0.0;
    i += 1;
    # Y axis
    dof_map.map[i] = (back_to_front_map[1][k]*3 - 1, back_to_front_map[2][k]*3 - 1);
    dof_map.offsets[i] = 0.0;
    i += 1;
    # Z axis
    dof_map.map[i] = (back_to_front_map[1][k]*3, back_to_front_map[2][k]*3);
    dof_map.offsets[i] = 0.0;
    i += 1;
end
i
dof_map.offsets

mpc = PenaltyMPC(dof_map, size(stiffness_matrix, 1), E * 1e4)
# rank(stiffness_matrix)
stiffness_matrix += mpc.penalty_matrix * mpc.penalty
rhs_vector += mpc.penalty_forces * mpc.penalty

##
# rank(stiffness_matrix)
central_node = get_nodes_by_tag(model, "center")
for i in -2:0
    stiffness_matrix[central_node[1]*3 + i, central_node[1]*3 + i] += mpc.penalty
end


# rank(stiffness_matrix)
# rhs_vector
using IterativeSolvers
# prob = LinearProblem(stiffness_matrix, rhs_vector)

# linsolve = init(prob);
# sol1 = solve(linsolve)

du = cg(stiffness_matrix, rhs_vector)
# rhs_vector
##
res_file = "results_new" |> results
uh_lin = FEFunction(trial_space_U,du)
writevtk(Ω, res_file ,  cellfields=["uh"=>uh_lin, "stress"=>σ∘ε(uh_lin)])


##


function get_nodes_by_tag(model::DiscreteModel, tag::String)::Vector{Int64}
    labels = get_face_labeling(model)

    dim = 0 # we will look for nodes which has dimension of 0
    
    dface_to_entity = get_face_entity(labels, dim)    
    dface_to_isontag = BitVector(undef, num_faces(labels,dim))
    tag_entities = get_tag_entities(labels, tag)

    for entity in tag_entities
        for i in eachindex(dface_to_entity)
            dface_to_isontag[i] =  dface_to_entity[i] == entity
        end
    end
    return findall(dface_to_isontag)
end

labels = get_face_labeling(model)
dface_to_entity = get_face_entity(labels, 0)


tag_entities = get_tag_entities(labels, "left_face")