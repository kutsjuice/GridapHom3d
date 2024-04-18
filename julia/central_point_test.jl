using Gmsh

gmsh.initialize()

gmsh.model.add("t18")
model = gmsh.model;

model.occ.addPoint(0, 0, 0, 0.1, 100);
model.occ.addPoint(0, 1, 0, 0.1, 101);
model.occ.addPoint(1, 1, 0, 0.1, 102);
model.occ.addPoint(1, 0, 0, 0.1, 103);
model.occ.addPoint(0.5, 0.5, 0, 0.01, 104);
model.occ.synchronize()


model.occ.addLine(100, 101, 200)
model.occ.addLine(101, 102, 201)
model.occ.addLine(102, 103, 202)
model.occ.addLine(103, 100, 203)


model.occ.synchronize()

model.occ.addCurveLoop([200, 201, 202, 203], 210)

model.occ.addPlaneSurface([210], 300)

model.occ.synchronize()

model.mesh.embed(0, [104], 2, 300 )
model.mesh.generate(2)
if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end

gmsh.finalize()