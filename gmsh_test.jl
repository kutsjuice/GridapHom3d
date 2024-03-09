# ------------------------------------------------------------------------------
#
#  Gmsh Julia tutorial 18
#
#  Periodic meshes
#
# ------------------------------------------------------------------------------

# Periodic meshing constraints can be imposed on surfaces and curves.

using Gmsh


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

next(ind, size) = ind==size ? 1 : ind+1;

function add_normal_periodic_cylinder(model, dir, pos, length, diam)
    
    cyl_sp = Vector{Float64}(undef, 3); # cylinder start point
    cyl_ep = Vector{Float64}(undef, 3); # cylinder end point

    cyl_sp[dir] = -length/2;
    cyl_ep[dir] = length;
    
    d2 = next(dir, 3);
    cyl_sp[d2]  = pos[1];
    cyl_ep[d2] = 0;

    d3 = next(d2, 3);
    cyl_sp[d3] =  pos[2];
    cyl_ep[d3] = 0;

    cyl = model.occ.addCylinder(cyl_sp..., cyl_ep..., diam);

    f_bbox = create_plate_selection_bbox(dir, -length/2, diam*2);

    s_bbox = create_plate_selection_bbox(dir, length/2, diam*2);
    f_bbox[[d2, d2+3]] .+= pos[1];
    s_bbox[[d2, d2+3]] .+= pos[1];
    f_bbox[[d3, d3+3]] .+= pos[2];
    s_bbox[[d3, d3+3]] .+= pos[2];
    model.occ.synchronize();
    first = model.getEntitiesInBoundingBox(f_bbox..., 2);
    second = model.getEntitiesInBoundingBox(s_bbox..., 2);
    println(model.getBoundingBox(3, cyl))
    println(f_bbox)
    println(s_bbox)
    println(first);
    println(second);
    first_tag = first[1][2];
    second_tag = second[1][2];

    translation = Vector{Float64}([
    1, 0, 0, 0, 
    0, 1, 0, 0, 
    0, 0, 1, 0, 
    0, 0, 0, 1
    ]);
    translation[dir*4] = length;

    # model.mesh.setPeriodic(2,   [second_tag], [first_tag], translation);
    
    return cyl;
end

function make_periodic_cell(model, cell_size)
    for d1 in [1,2,3]
        # d2 = next(d1, 3);
        # d3 = next(d3, 3);

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
end

gmsh.initialize()

gmsh.model.add("t18")
model = gmsh.model;
# Let's use the OpenCASCADE geometry kernel to build two geometries.

# The first geometry is very simple: a unit cube with a non-uniform mesh size
# constraint (set on purpose to be able to verify visually that the periodicity
# constraint works!):

diam = 0.1;

cell_size = 0.5;


add_normal_periodic_cylinder(model, 1, [-0.05, 0.05], cell_size, diam);
add_normal_periodic_cylinder(model, 2, [-0.05, 0.05], cell_size, diam);
add_normal_periodic_cylinder(model, 3, [-0.05, 0.05], cell_size, diam);

model.occ.synchronize()

dimtags = model.getEntities(3)
model.occ.fuse(dimtags[1:end-1], dimtags[end:end])
model.occ.synchronize()

make_periodic_cell(model, cell_size);
model.occ.synchronize()

MSFC = 3;
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", MSFC);
gmsh.option.setNumber("Mesh.MeshSizeMax", diam/5);
gmsh.option.setNumber("Mesh.MeshSizeMin", diam/10);

model.mesh.generate(3)
gmsh.write("periodic.vtk")
gmsh.write("periodic.msh")

# Launch the GUI to see the results:
# if !("-nopopup" in ARGS)
#     gmsh.fltk.run()
# end

gmsh.finalize()
