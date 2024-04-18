using Gridap
using Gridap.Fields
using Gridap.CellData
using Gridap.Geometry
using Gridap.ReferenceFEs
using FillArrays

pmin = Point(0,0,0); pmax = Point(1, 1, 1)

n = 10; cells = (n, n, n)
model = CartesianDiscreteModel(pmin, pmax, cells)

Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)

q = Point{3, Float64}[(.5, .5, .5)]
s = Point{2, Float64}[(.5, .5)]


cell_q = Fill(q, num_cells(Ω))
cell_s = Fill(s, num_cells(Γ))

x_Ω = CellPoint(cell_q, Ω, ReferenceDomain())
x_Γ = CellPoint(cell_s, Γ, ReferenceDomain())

ffun(x) = x[1] + x[2] + x[3]

f = GenericField(ffun)

cell_f = Fill(f, num_cells(Ω))
f_Ω = GenericCellField(cell_f, Ω, PhysicalDomain())


fx_Ω = f_Ω(x_Ω)

grid = get_grid(Ω)
conn = get_grid_topology(model)

get_node_coordinates(Ω)[1]