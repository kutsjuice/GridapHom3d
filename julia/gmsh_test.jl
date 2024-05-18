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

gmsh.model.add("three_cross")
model = gmsh.model;

diam = 0.1;

cell_size = 0.5;

offset = 0.05
add_normal_periodic_cylinder(model, 1, [-offset, offset], cell_size, diam);
add_normal_periodic_cylinder(model, 2, [-offset, offset], cell_size, diam);
add_normal_periodic_cylinder(model, 3, [-offset, offset], cell_size, diam);

model.occ.synchronize()

dimtags = model.getEntities(3)
model.occ.fuse(dimtags[1:end-1], dimtags[end:end])
model.occ.synchronize()

EPS = 1e-2
min_value = -diam/2 - EPS - offset
max_value = diam/2 + EPS + offset
center_intesection_bbox = [
    min_value, 
    min_value, 
    min_value, 
    max_value, 
    max_value, 
    max_value
]

model.getEntitiesInBoundingBox(center_intesection_bbox..., 2)


# Launch the GUI to see the results:
if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end


make_periodic_cell(model, cell_size);
model.occ.synchronize()

volume = model.getEntities(3)

model.addPhysicalGroup(volume[1][1], [volume[1][2]], 300, "cell")

central_point = model.occ.addPoint(0, 0, 0, 0.0001, 100);
model.occ.synchronize()

model.addPhysicalGroup(0, [100], -1, "center")

model.mesh.embed(0, [100], volume[1][1], volume[1][2])


# the rule is as follow:
#
#  -   LEFT = [X-]
#  -  RIGHT = [X+]
#
#  - BOTTOM = [Y-]
#  -    TOP = [Y+]
#
#  -   BACK = [Z-]
#  -  FRONT = [Z+]

# select left and right surfaces along X axis
left_face_bbox = [
    -cell_size/2 - 1e-3,
    -cell_size/2 - 1e-3,
    -cell_size/2 - 1e-3,

    -cell_size/2 + 1e-3, # <<==
    cell_size/2 + 1e-3,
    cell_size/2 + 1e-3,
]

right_face_bbox = [
    cell_size/2 - 1e-3, # <<==
    -cell_size/2 - 1e-3,
    -cell_size/2 - 1e-3,

    cell_size/2 + 1e-3,
    cell_size/2 + 1e-3,
    cell_size/2 + 1e-3,
]

# select bottom and top surfaces along Y axis
bottom_face_bbox = [
    -cell_size/2 - 1e-3,
    -cell_size/2 - 1e-3,
    -cell_size/2 - 1e-3,

    cell_size/2 + 1e-3,  
    -cell_size/2 + 1e-3,  # <<==
    cell_size/2 + 1e-3,
]

top_face_bbox = [
    -cell_size/2 - 1e-3,
    cell_size/2 - 1e-3, # <<==
    -cell_size/2 - 1e-3,

    cell_size/2 + 1e-3,
    cell_size/2 + 1e-3,
    cell_size/2 + 1e-3,
]

# select front and back surfaces along Z axis
back_face_bbox = [
    -cell_size/2 - 1e-3,
    -cell_size/2 - 1e-3,
    -cell_size/2 - 1e-3,

    cell_size/2 + 1e-3,  
    cell_size/2 + 1e-3,
    -cell_size/2 + 1e-3,  # <<==
]

front_face_bbox = [
    -cell_size/2 - 1e-3,
    -cell_size/2 - 1e-3,
    cell_size/2 - 1e-3, # <<==

    cell_size/2 + 1e-3,
    cell_size/2 + 1e-3,
    cell_size/2 + 1e-3,
]



left_face   = model.getEntitiesInBoundingBox(left_face_bbox..., 2);
right_face  = model.getEntitiesInBoundingBox(right_face_bbox..., 2);

bottom_face = model.getEntitiesInBoundingBox(bottom_face_bbox..., 2);
top_face    = model.getEntitiesInBoundingBox(top_face_bbox..., 2);

back_face   = model.getEntitiesInBoundingBox(back_face_bbox..., 2)
front_face  = model.getEntitiesInBoundingBox(front_face_bbox..., 2)

model.addPhysicalGroup(left_face[1][1], [left_face[1][2]], -1, "left_face")
model.addPhysicalGroup(right_face[1][1], [right_face[1][2]], -1, "right_face")

model.addPhysicalGroup(bottom_face[1][1], [bottom_face[1][2]], -1, "bottom_face")
model.addPhysicalGroup(top_face[1][1], [top_face[1][2]], -1, "top_face")

model.addPhysicalGroup(back_face[1][1], [back_face[1][2] ], -1, "back_face" )
model.addPhysicalGroup(front_face[1][1], [front_face[1][2]], -1, "front_face")




left_face   = model.getEntitiesInBoundingBox(left_face_bbox..., 1)
right_face  = model.getEntitiesInBoundingBox(right_face_bbox..., 1)

bottom_face = model.getEntitiesInBoundingBox(bottom_face_bbox..., 1)
top_face    = model.getEntitiesInBoundingBox(top_face_bbox..., 1)

back_face   = model.getEntitiesInBoundingBox(back_face_bbox..., 1)
front_face  = model.getEntitiesInBoundingBox(front_face_bbox..., 1)

model.addPhysicalGroup(left_face[1][1], [left_face[1][2]], -1, "left_face")
model.addPhysicalGroup(right_face[1][1], [right_face[1][2]], -1, "right_face")

model.addPhysicalGroup(bottom_face[1][1], [bottom_face[1][2]], -1, "bottom_face")
model.addPhysicalGroup(top_face[1][1], [top_face[1][2]], -1, "top_face")

model.addPhysicalGroup(back_face[1][1], [back_face[1][2] ], -1, "back_face" )
model.addPhysicalGroup(front_face[1][1], [front_face[1][2]], -1, "front_face")


left_face   = model.getEntitiesInBoundingBox(left_face_bbox..., 0)
right_face  = model.getEntitiesInBoundingBox(right_face_bbox..., 0)

bottom_face = model.getEntitiesInBoundingBox(bottom_face_bbox..., 0)
top_face    = model.getEntitiesInBoundingBox(top_face_bbox..., 0)

back_face   = model.getEntitiesInBoundingBox(back_face_bbox..., 0)
front_face  = model.getEntitiesInBoundingBox(front_face_bbox..., 0)

model.addPhysicalGroup(left_face[1][1], [left_face[1][2]], -1, "left_face")
model.addPhysicalGroup(right_face[1][1], [right_face[1][2]], -1, "right_face")

model.addPhysicalGroup(bottom_face[1][1], [bottom_face[1][2]], -1, "bottom_face")
model.addPhysicalGroup(top_face[1][1], [top_face[1][2]], -1, "top_face")

model.addPhysicalGroup(back_face[1][1], [back_face[1][2] ], -1, "back_face" )
model.addPhysicalGroup(front_face[1][1], [front_face[1][2]], -1, "front_face")

# model.addPhysicalGroup(0, [central_point])

MSFC = 3;
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", MSFC);
gmsh.option.setNumber("Mesh.MeshSizeMax", diam/5);
gmsh.option.setNumber("Mesh.MeshSizeMin", diam/10);
gmsh.option.setNumber("Mesh.ElementOrder", 2)#колличество элементов на 2pi радиан



node = model.getEntitiesInBoundingBox(-1e-3, -1e-3, -1e-3, 1e-3, 1e-3, 1e-3) 

model.mesh.generate(3)
input(filename::String) = "geometry/" * filename;


# p_rf = model.mesh.get_periodic_nodes(right_face[1][1], right_face[1][2])
# p_lf = model.mesh.get_periodic_nodes(left_face[1][1], left_face[1][2])

# p_rf = model.mesh.get_periodic_nodes(right_face[1][1], right_face[1][2])
# p_lf = model.mesh.get_periodic_nodes(left_face[1][1], left_face[1][2])

# p_rf = model.mesh.get_periodic_nodes(right_face[1][1], right_face[1][2])
# p_lf = model.mesh.get_periodic_nodes(left_face[1][1], left_face[1][2])

# print(length(p_rf[2]))

gmsh.write("periodic.vtk" |> input)
gmsh.write("periodic.msh" |> input)


# using GridapGmsh



# Launch the GUI to see the results:
# if !("-nopopup" in ARGS)
#     gmsh.fltk.run()
# end

gmsh.finalize()

# ppp = p_rf[2] |> (x) -> map(Int, x)
# vvv = p_rf[3] |> (x) -> map(Int, x)



# sort(ppp)
# sort(vvv)