# ------------------------------------------------------------------------------
#
#  Gmsh Julia tutorial 18
#
#  Periodic meshes
#
# ------------------------------------------------------------------------------

# Periodic meshing constraints can be imposed on surfaces and curves.

using Gmsh

gmsh.initialize()

gmsh.model.add("t18")

# Let's use the OpenCASCADE geometry kernel to build two geometries.

# The first geometry is very simple: a unit cube with a non-uniform mesh size
# constraint (set on purpose to be able to verify visually that the periodicity
# constraint works!):

diam = 0.1;

cell_size = 0.5;

fabric = gmsh.model.occ;

fabric.addCylinder(-cell_size/2, -diam, 0, cell_size, 0,0,diam/2);
fabric.addCylinder(-cell_size/2, diam,0, cell_size, 0,0,diam/2);
fabric.synchronize();

fabric.getEntities(2)

eps = 0.01;

#select left boundaries
lbb = [-cell_size/2 - eps, -cell_size, -cell_size, -cell_size/2 + eps, cell_size, cell_size]
left = gmsh.model.getEntitiesInBoundingBox(lbb...,2)
#select right boundaries
rbb = [cell_size/2 - eps, -cell_size, -cell_size, cell_size/2 + eps, cell_size, cell_size]
right = gmsh.model.getEntitiesInBoundingBox(rbb..., 2)

translation = [
    1, 0, 0, cell_size, 
    0, 1, 0, 0, 
    0, 0, 1, 0, 
    0, 0, 0, 1
    ];

left_faces = [dimtag[2] for dimtag in left]
right_faces = [dimtag[2] for dimtag in right]


gmsh.model.mesh.setPeriodic(2,  right_faces, left_faces, translation)



fabric.synchronize()

MSFC = 3;
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", MSFC);
gmsh.option.setNumber("Mesh.MeshSizeMax", diam/5);
gmsh.option.setNumber("Mesh.MeshSizeMin", diam/10);

gmsh.model.mesh.generate(3)
gmsh.write("periodic.vtk")

# Launch the GUI to see the results:
# if !("-nopopup" in ARGS)
#     gmsh.fltk.run()
# end

gmsh.finalize()
