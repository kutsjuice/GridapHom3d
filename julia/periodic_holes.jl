using Gmsh

input(filename::String) = "geometry/" * filename;
results(filename::String) = "results/" * filename;


include("pbc_gmsh_utils.jl")
gmsh.finalize()
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0);

gmsh.model.add("t18")
model = gmsh.model;
m = 1;
mm = 1e-3;
μm = 1e-6;
nm = 1e-9;
cell_size = 1m;
φ = 0.6
dₚ = 0.2*cell_size
V_RVE = cell_size^3
Vₚ = 4/3*π*dₚ^3
N_approx = 100;ceil(Int, V_RVE / Vₚ)
order = 1

box = model.occ.add_box(- cell_size/2, - cell_size/2, - cell_size/2, cell_size, cell_size, cell_size, -1)

model.occ.synchronize()
# xc, yc, zc = rand(3) * cell_size .- cell_size/2
# cs = cell_size

# add_periodically(-0.5, 0., 0., 0.1, 1)


pores = Vector{Tuple{Int64, Int64}}()

for i in 1:N_approx
    pore_center = rand(3) * cell_size .- cell_size/2
    pores_buf = add_periodical_sphere(pore_center..., dₚ/2, cell_size);
    push!(pores, pores_buf...)
end

pores
# central_point = model.occ.addPoint(0, 0, 0, 0.0001, 100);

model.occ.synchronize()
##

model.occ.cut([3,box], pores)
model.occ.synchronize()

# for d1 in 1:3
    d1 = 1
    f_bbox = create_plate_selection_bbox(d1, -cell_size/2, cell_size*2);
    s_bbox = create_plate_selection_bbox(d1, cell_size/2, cell_size*2);

    dim = 2
        first_dimtags = model.occ.getEntitiesInBoundingBox(f_bbox..., dim)
        second_dimtags = model.occ.getEntitiesInBoundingBox(s_bbox..., dim)

        first_tags = [dimtag[2] for dimtag in first_dimtags];
        second_tags = [dimtag[2] for dimtag in second_dimtags];

        tr_mat = Vector{Float64}([
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
        ]);
        tr_mat[d1*4] = cell_size;

        model.mesh.setPeriodic(dim, second_tags, first_tags, tr_mat);
    

# end
EPS = cell_size/100
min_lim = - cell_size/2 - EPS
max_lim = cell_size/2 + EPS
##
model.occ.synchronize()

body = model.occ.getEntities(3)
model.occ.getBoundingBox(body[1]...)
##
vₛ = gmsh.model.occ.getMass(body[1]...)

gmsh.option.setNumber("Mesh.MeshSizeMax", cell_size/20);
gmsh.option.setNumber("Mesh.MeshSizeMin", dₚ/10);
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 20)
gmsh.option.setNumber("Mesh.ElementOrder", order)
model.mesh.generate(3)
if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end
msh_file = "pores.msh" |> input
gmsh.write("pores.vtk" |> input)
gmsh.write(msh_file)

gmsh.finalize()
##
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
