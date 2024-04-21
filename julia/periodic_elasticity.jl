using Gridap
using GridapGmsh
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.Visualization
using SparseArrays

using LinearSolve
# using LinearAlgebra
# using LinearSolve

using Gmsh

input(filename::String) = "geometry/" * filename;
results(filename::String) = "results/" * filename;

include("penalty_mpc.jl")
include("pbc_mpc_utils.jl")

order = 2


function build_relation_map(first_nodes::Vector, second_nodes::Vector, model::DiscreteModel, relation::Function)
    
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

function pointhash(point)
    return point[1] + point[2] * 1e3 + point[3] * 1e6 
end

function build_relation_map_by_hash(first_nodes::Vector, second_nodes::Vector, model::DiscreteModel, relation::Function, hash::Function)
    
    @assert length(first_nodes) == length(second_nodes) "Size of node sets doesn't match"

    node_coordinates = get_node_coordinates(model);

    first_hashes = lazy_map(hash, node_coordinates[first_nodes]);
    first_nodes_sorted = sortperm(first_hashes);
    second_hashes = lazy_map(hash, node_coordinates[second_nodes]);
    second_nodes_sorted = sortperm(second_hashes);

    check = BitVector(undef, length(first_nodes))
    # ch = Vector{Float64}(undef, length(first_nodes))
    for (i, (first_ind, second_ind)) in enumerate(zip(first_nodes[first_nodes_sorted], second_nodes[second_nodes_sorted]))
        check[i] = relation(node_coordinates[first_ind], node_coordinates[second_ind]) < 1e-5
    end

    # return ch
    
    @assert all(check) "Node sets don't match the relationship"

    return (first_nodes[first_nodes_sorted], second_nodes[second_nodes_sorted])
end


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


gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0);

gmsh.model.add("t18")
model = gmsh.model;

cell_size = 0.5;

model.occ.add_box(- cell_size/2, - cell_size/2, - cell_size/2, cell_size, cell_size, cell_size, -1)

model.occ.synchronize()

for d1 in 1:3
    f_bbox = create_plate_selection_bbox(d1, -cell_size/2, cell_size*2);
    s_bbox = create_plate_selection_bbox(d1, cell_size/2, cell_size*2);

    first_dimtags = model.occ.getEntitiesInBoundingBox(f_bbox..., 2);
    second_dimtags = model.occ.getEntitiesInBoundingBox(s_bbox..., 2);

    first_tags = [dimtag[2] for dimtag in first_dimtags];
    second_tags = [dimtag[2] for dimtag in second_dimtags];

    tr_mat = Vector{Float64}([
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    ]);
    tr_mat[d1*4] = cell_size;

    model.mesh.setPeriodic(2, second_tags, first_tags, tr_mat);
end

central_point = model.occ.addPoint(0, 0, 0, 0.0001, 100);
model.occ.synchronize()

model.addPhysicalGroup(0, [100], -1, "center")

volume = model.getEntities(3)
model.addPhysicalGroup(volume[1][1], [volume[1][2]], 300, "cell")


model.mesh.embed(0, [100], volume[1][1], volume[1][2])


gmsh.option.setNumber("Mesh.MeshSizeMax", 0.03);
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.01);
gmsh.option.setNumber("Mesh.ElementOrder", order)


model.mesh.generate(3)
if !in("geometry", readdir())
    mkdir("geometry")
end
msh_file = "box.msh" |> input
gmsh.write("box.vtk" |> input)
gmsh.write(msh_file)
gmsh.finalize()

global fe_model 

for i in 1:2
    try
        fe_model = GmshDiscreteModel(msh_file, renumber=false)
    catch Error
    end
end

vtk_file = "imported_mesh" |> input;
writevtk(fe_model, vtk_file);

x_neg_bbox = create_plate_selection_bbox(1, -cell_size/2, cell_size*2);
x_pos_bbox = create_plate_selection_bbox(1, cell_size/2, cell_size*2);

y_neg_bbox = create_plate_selection_bbox(2, -cell_size/2, cell_size*2)
y_pos_bbox = create_plate_selection_bbox(2, cell_size/2, cell_size*2);

z_neg_bbox = create_plate_selection_bbox(3, -cell_size/2, cell_size*2);
z_pos_bbox = create_plate_selection_bbox(3, cell_size/2, cell_size*2);


x_neg_bound_nodes = get_nodes_in_bounding_box(fe_model, x_neg_bbox);
x_pos_bound_nodes = get_nodes_in_bounding_box(fe_model, x_pos_bbox);

y_neg_bound_nodes = get_nodes_in_bounding_box(fe_model, y_neg_bbox)
y_pos_bound_nodes = get_nodes_in_bounding_box(fe_model, y_pos_bbox)

z_neg_bound_nodes = get_nodes_in_bounding_box(fe_model, z_neg_bbox);
z_pos_bound_nodes = get_nodes_in_bounding_box(fe_model, z_pos_bbox);


const E = 70.0e9
const ν = 0.33
const λ = (E*ν)/((1+ν)*(1-2*ν))
const μ = E/(2*(1+ν))


function stress(x, strain)
    elastic_module = 70.0e9 * (1 + 0.2 * sin(x[1] * 2π / cell_size) * sin(x[2] * 2π / cell_size) * sin(x[3] * 2π / cell_size)) ;
    poisson_coeff  = 0.3 * (1 + 0.4 * sin(x[1] * 2π / cell_size) * sin(x[2] * 2π / cell_size));
    lame_paramter1 = (elastic_module * poisson_coeff) / ( (1 + poisson_coeff) * (1 - 2*poisson_coeff) );
    lame_paramter2 = elastic_module / ( 2 * (1 + poisson_coeff) );

    return lame_paramter1 * tr(strain) * one(strain) + 2 * lame_paramter2 * strain;
end

# σ(ε) = λ*tr(ε)*one(ε) + 2*μ*ε

degree = order*2
Ω = Triangulation(fe_model)

dΩ = Measure(Ω,degree)

reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
test_space_V = TestFESpace(fe_model,reffe,conformity=:H1)#,dirichlet_tags = ["center"])
trial_space_U = TrialFESpace(test_space_V)

# weak form 
px = get_physical_coordinate(Ω)

bilinear_form_stiffness(u,v) = ∫( ε(v) ⊙ (stress∘(px, ε(u))) )*dΩ
linear_form_stiffness(v) = 0; #∫(dot(f, v))*dΩ

operator_K = AffineFEOperator(bilinear_form_stiffness, linear_form_stiffness, trial_space_U, test_space_V)

stiffness_matrix = get_matrix(operator_K)
rhs_vector = get_vector(operator_K)

function x_periodic_relation(negative::VectorValue, positive::VectorValue)
    return norm( positive - negative - VectorValue(0.5, 0.0, 0.0))
end

function y_periodic_relation(negative::VectorValue, positive::VectorValue)
    return norm(positive - negative - VectorValue(0.0, 0.5, 0.0))
end

function z_periodic_relation(negative::VectorValue, positive::VectorValue)
    return norm(positive - negative - VectorValue(0.0, 0.0, 0.5))
end

x_neg_to_pos_map = build_relation_map_by_hash(x_neg_bound_nodes, x_pos_bound_nodes, fe_model, x_periodic_relation, pointhash)
y_neg_to_pos_map = build_relation_map_by_hash(y_neg_bound_nodes, y_pos_bound_nodes, fe_model, y_periodic_relation, pointhash)
z_neg_to_pos_map = build_relation_map_by_hash(z_neg_bound_nodes, z_pos_bound_nodes, fe_model, z_periodic_relation, pointhash)

n1 = length(x_neg_to_pos_map[1])
n2 = length(y_neg_to_pos_map[1])
n3 = length(z_neg_to_pos_map[1])

# get_nodes_by_tag(model, "left_face")

node_to_dofs = Vector{Tuple{Int64, Int64, Int64}}(undef, num_nodes(fe_model))
node_is_set = BitVector(undef, num_nodes(fe_model)); fill!(node_is_set, 0);
cell_dofs = get_cell_dof_ids(trial_space_U)
cell_nodes = get_cell_node_ids(fe_model)

for (dofs, nodes) in zip(cell_dofs, cell_nodes)
    for (i, node) in enumerate(nodes)
        if !node_is_set[node]
            node_to_dofs[node] = (dofs[i], dofs[i + length(nodes)], dofs[i + 2*length(nodes)]);
            node_is_set[node] = 1;
        else
            # println(dofs')
            # println(nodes')
            if (dofs[i], dofs[i + length(nodes)], dofs[i + 2*length(nodes)]) != node_to_dofs[node]
                println(dofs')
                println(nodes')
            end
        end
    end
end

i = 1;
dof_map = Dof2DofMap((n1+n2+n3)*3)

(n1+n2+n3)*3
for k in eachindex(x_neg_to_pos_map[1])
    # X axis
    dof_map.map[i] = (
        node_to_dofs[x_neg_to_pos_map[1][k]][1], 
        node_to_dofs[x_neg_to_pos_map[2][k]][1]
        );
    dof_map.offsets[i] = 0.0;
    i += 1;
    # Y axis
    dof_map.map[i] = (
        node_to_dofs[x_neg_to_pos_map[1][k]][2], 
        node_to_dofs[x_neg_to_pos_map[2][k]][2]
        );
    dof_map.offsets[i] = 0.05;
    i += 1;
    # Z axis
    dof_map.map[i] = (
        node_to_dofs[x_neg_to_pos_map[1][k]][3], 
        node_to_dofs[x_neg_to_pos_map[2][k]][3]
        );
    dof_map.offsets[i] = 0.0;
    i += 1;
end
for k in eachindex(y_neg_to_pos_map[1])
    # X axis
    dof_map.map[i] = (
        node_to_dofs[y_neg_to_pos_map[1][k]][1], 
        node_to_dofs[y_neg_to_pos_map[2][k]][1]
        );
    dof_map.offsets[i] = 0.05;
    i += 1;
    # Y axis
    dof_map.map[i] = (
        node_to_dofs[y_neg_to_pos_map[1][k]][2], 
        node_to_dofs[y_neg_to_pos_map[2][k]][2]
        );
    dof_map.offsets[i] = 0.0;
    i += 1;
    # Z axis
    dof_map.map[i] = (
        node_to_dofs[y_neg_to_pos_map[1][k]][3], 
        node_to_dofs[y_neg_to_pos_map[2][k]][3]
        );
    dof_map.offsets[i] = 0.0;
    i += 1;
end
for k in eachindex(z_neg_to_pos_map[1])
    # X axis
    dof_map.map[i] = (
        node_to_dofs[z_neg_to_pos_map[1][k]][1], 
        node_to_dofs[z_neg_to_pos_map[2][k]][1]
        );
    dof_map.offsets[i] = 0.0;
    i += 1;
    # Y axis
    dof_map.map[i] = (
        node_to_dofs[z_neg_to_pos_map[1][k]][2], 
        node_to_dofs[z_neg_to_pos_map[2][k]][2]
        );
    dof_map.offsets[i] = 0.0;
    i += 1;
    # Z axis
    dof_map.map[i] = (
        node_to_dofs[z_neg_to_pos_map[1][k]][3], 
        node_to_dofs[z_neg_to_pos_map[2][k]][3]
        );
    dof_map.offsets[i] = 0.0;
    i += 1;
end
i
dof_map.offsets
# stiffness_matrix
dof_map

mpc = PenaltyMPC(dof_map, size(stiffness_matrix, 1), E * 1e4)
# rank(stiffness_matrix)
stiffness_matrix += mpc.penalty_matrix * mpc.penalty
rhs_vector += mpc.penalty_forces * mpc.penalty


# rank(stiffness_matrix)
central_node = get_nodes_by_tag(fe_model, "center")
for i in -2:0
    stiffness_matrix[central_node[1]*3 + i, central_node[1]*3 + i] += mpc.penalty
end

#
# rank(stiffness_matrix)
# rhs_vector
using IterativeSolvers
using LinearSolve
using AlgebraicMultigrid
using GeometricMultigrid
# import IterativeSolvers: cg

ml = ruge_stuben(stiffness_matrix) 
p = aspreconditioner(ml)
# du = cg(stiffness_matrix, rhs_vector, Pl = p)

# du = solve(stiffness_matrix, rhs_vector, RugeStubenAMG(), maxiter = 1, abstol = 1e-6)
# du, it = mg(stiffness_matrix, rhs_vector)
prob = LinearProblem(stiffness_matrix, rhs_vector)
sol = solve(prob, KrylovJL_CG(), Pl = p)
# sol.u
# sol1 = solve(prob)
du = sol.u
#
# du = cg(stiffness_matrix, rhs_vector)
# rhs_vector
#
res_file = "cube" |> results
uh_lin = FEFunction(trial_space_U,du)
writevtk(Ω, res_file ,  cellfields=["uh"=>uh_lin, "stress"=>stress∘(px, ε(uh_lin))])


##
# cell_node_ids = get_cell_node_ids(Ω)
# cell_dof_ids = get_cell_dof_ids(test_space_V)
# cell_node_ids[1]'
# cell_dof_ids[1]'?