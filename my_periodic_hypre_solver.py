import numpy as np;
import ufl;
from inspect import currentframe, getframeinfo

cf = currentframe()
filename = getframeinfo(cf).filename

from typing import Union
from pathlib import Path

from mpi4py import MPI;
from petsc4py import PETSc

from dolfinx import mesh, fem, plot, io, la;

from dolfinx.fem.petsc import apply_lifting, assemble_matrix, assemble_vector, set_bc
import dolfinx_mpc.utils
from dolfinx_mpc import LinearProblem, MultiPointConstraint
from dolfinx.io import XDMFFile, gmshio;
import dolfinx.geometry as geo;
import gmsh;

import sys;

print("RUNNING: ")#, "in/"+sys.argv[1])
D_TYPE = PETSc.ScalarType

ν = 0.43;
E = 1.8e9;

λ = ν*E/(1+ν)/(1-2*ν);
μ = E/2/(1+ν);

def epsilon(u):
    return ufl.sym(ufl.grad(u)); # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)
def sigma(u):
    return λ * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2*μ*epsilon(u);


def build_nullspace(V: fem.VectorFunctionSpace):
    """Build  PETSc nullspace for 3D elasticity"""
    
    # Create vector that span the nullspace
    bs = V.dofmap.index_map_bs;
    length0 = V.dofmap.index_map.size_local;
    length1 = length0 + V.dofmap.index_map.num_ghosts;
    basis = [np.zeros(bs * length1, dtype = D_TYPE) for i in range(6)];
    
    # Get dof indices for each subspace (x, y and z dofs)
    dofs = [V.sub(i).dofmap.list.flatten() for i in range(3)];
    
    # Set the three translational rigid body modes
    for i in range(3):
        basis[i][dofs[i]] = 1.0;
    
    # Set the three rotational rigid body modes
    x = V.tabulate_dof_coordinates();
    dofs_block = V.dofmap.list.flatten();
    x0, x1, x2 = x[dofs_block, 0], x[dofs_block, 1], x[dofs_block, 2];
    
    basis[3][dofs[0]] = -x1;
    basis[3][dofs[1]] = x0;
    basis[4][dofs[0]] = x2;
    basis[4][dofs[2]] = -x0;
    basis[5][dofs[2]] = x1;
    basis[5][dofs[1]] = -x2;
    
    # Create PETSc Vec objects (excluding ghosts) and normalise
    basis_petsc = [PETSc.Vec().createWithArray(x[:bs*length0], bsize=3, comm=V.mesh.comm) for x in basis]
    la.orthonormalize(basis_petsc);
    assert la.is_orthonormal(basis_petsc);
    
    #Create and return a PETSc nullspace
    return PETSc.NullSpace().create(vectors=basis_petsc);


def main():
    
    ORDER = 1

    
    ## Setting up gmsh properties
    gmsh.initialize()

    # Choose if Gmsh output is verbose
    gmsh.option.setNumber("General.Terminal", 0)

    # Set elements order to the specified one
    gmsh.option.setNumber("Mesh.ElementOrder", ORDER)
    # Set elements size
    # gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 5) # uncomment to use for mesh refinement dependending from its surface curvature
    gmsh.option.setNumber("Mesh.MeshSizeMax", 5e-2)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 1e-2)

    # Set threads number for distrebuted meshing
    # gmsh.option.setNumber("Mesh.MaxNumThreads3D", 4)

    # Set mesh algorithm (default is Delaunay triangulation)
    # see https://gmsh.info/doc/texinfo/gmsh.html#Choosing-the-right-unstructured-algorithm
    gmsh.option.setNumber("Mesh.Algorithm3D", 3)

    # gmsh.option.setNumber("Mesh.RecombinationAlgorithm",3)
    # gmsh.option.setNumber("Mesh.Recombine3DAll",1)

    # Set the usage of hexahedron elements 
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 0)
    ## Importing RVE geometry

    # gmsh.open("in/"+sys.argv[1]);
    gmsh.open("periodic.msh");
    
    model = gmsh.model()
    # model.add("main_domain")
    model_name = model.getCurrent()
    tags = [dimtag[1] for dimtag in model.get_entities(3)]

    model.add_physical_group(dim=3, tags=tags)


    # Synchronize OpenCascade representation with gmsh model
    model.occ.synchronize()

    # Generate the mesh
    # model.mesh.generate(2)
    # model.mesh.recombine()
    model.mesh.generate(dim=3)

    bbox = [np.Inf,
            np.Inf,
            np.Inf,
            -np.Inf,
            -np.Inf,
            -np.Inf]
    for tag in tags:
        buf_bbox = model.get_bounding_box(3, tag)
        for i in range(3):
            if bbox[i] > buf_bbox[i]:
                bbox[i] = buf_bbox[i]
        for j in range(3,6):
            if bbox[j] < buf_bbox[j]:
                bbox[j] = buf_bbox[j]


    # Create a DOLFINx mesh (same mesh on each rank)
    # print("MPI.COMM_SELF====",MPI.COMM_SELF)
    msh, cell_markers, facet_markers = gmshio.model_to_mesh(model, MPI.COMM_WORLD,0)
    # msh, cell_markers, facet_markers = gmshio.read_from_msh("in/"+sys.argv[1], MPI.COMM_WORLD, 0, gdim=3)
    msh.name = "Box"
    cell_markers.name = f"{msh.name}_cells"
    facet_markers.name = f"{msh.name}_facets"
    
    # Finalize gmsh to be able to use it again
    gmsh.finalize()
    print("MESH IMPORTED")
    print("NUMBER OF NODES:", msh.geometry.x.shape[0])
   
    eps = np.linalg.norm(np.array(bbox[0:3]) + np.array(bbox[3:]));
    
    unit_disp =np.mean(np.array(bbox[3:]) - np.array(bbox[:3]));
    
    fdim = msh.topology.dim - 1

    def left(x):
        return np.isclose(x[0], bbox[0], atol = eps);

    def right(x):
        return np.isclose(x[0], bbox[3], atol = eps);

    def bottom(x):
        return np.isclose(x[2], bbox[2], atol = eps);

    def top(x):
        return np.isclose(x[2], bbox[5], atol = eps);

    def front(x):
        return np.isclose(x[1], bbox[1], atol = eps);

    def back(x):
        return np.isclose(x[1], bbox[4], atol = eps);

    def KUBC(x, i, j, ud):
        values = np.zeros(x.shape);

        values[i,:] += 0.5*ud*(x[j])/(bbox[j+3] - bbox[j]);
        values[j,:] += 0.5*ud*(x[i])/(bbox[i+3] - bbox[i]);

        return values;

    # find all facets on top, bottom and left boundary
    left_facets = mesh.locate_entities_boundary(msh, fdim, left);
    right_facets = mesh.locate_entities_boundary(msh, fdim, right);
    bottom_facets = mesh.locate_entities_boundary(msh, fdim, bottom);
    top_facets = mesh.locate_entities_boundary(msh, fdim, top);
    front_facets = mesh.locate_entities_boundary(msh, fdim, front);
    back_facets = mesh.locate_entities_boundary(msh, fdim, back);
    
    marked_facets = np.hstack([left_facets, 
                               right_facets, 
                               bottom_facets,
                               top_facets,
                               front_facets,
                               back_facets,
                              ]);

    markers = np.hstack([np.full_like(left_facets, 1),
                         np.full_like(right_facets, 2),
                         np.full_like(bottom_facets, 3),
                         np.full_like(top_facets, 4),
                         np.full_like(front_facets, 5),
                         np.full_like(back_facets, 6),
                        ]);

    facets_order = np.argsort(marked_facets);

    facets_tags = mesh.meshtags(msh, 
                                fdim, 
                                marked_facets[facets_order],
                                markers[facets_order]);
    
    VF_space = fem.VectorFunctionSpace(msh, ("CG", ORDER))
    ds = ufl.Measure('ds', domain=msh, subdomain_data=facets_tags, metadata={'quadrature_degree': ORDER});    
    dx = ufl.Measure('dx', domain=msh, metadata={'quadrature_degree': ORDER});

    # Set displacements for Dirichlet BC

    ## You can choose to set the same displacement to all dofs or specific sisplacement to each dof.
    ## in the first case you can just specify ub_ as np.array([msh.topology.dim]) - 
    shifted_u = np.array([0.0, -0.1, 0.0])
    fixed_u = np.array([0.0, 0.0, 0.0])

    ## in the second case you need to provide np.array([x.shape]). For that you need to interpolate some function over the domain ass follow
    # ub_ = fem.Function(V);                        # <--- blanc displacements
    # full_bc = lambda x: KUBC(x, i,j, unit_disp);  # <--- functions which will be interpolated over domain
    # ub_.interpolate(full_bc);                     # <--- computed displacements

    shifted_dofs = fem.locate_dofs_topological(VF_space, 
                                               fdim,
                                               right_facets);
    fixed_dofs = fem.locate_dofs_topological(VF_space,
                                           fdim,
                                           left_facets);
    
    top_dofs = fem.locate_dofs_topological(VF_space, 
                                           fdim, 
                                           back_facets);
    bottom_dofs = fem.locate_dofs_topological(VF_space,
                                              fdim,
                                              front_facets);
    
    shifted_bc = fem.dirichletbc(shifted_u, shifted_dofs, VF_space);
    fixed_bc = fem.dirichletbc(fixed_u, fixed_dofs, VF_space);

    bc_combined = [fixed_bc, shifted_bc]

    cords_of_top_dofs = msh.geometry.x[top_dofs];
    cords_of_bottom_dofs = msh.geometry.x[bottom_dofs];

    for bottom

    ## let the back nodes are primary and the front are secondary
    def locate_secondary(x):
        return front(x)

    # bbox[1] - lower boundary
    # bbox[4] - upper boundary
    def periodic_relation(x):
        out_x = np.zeros_like(x)
        out_x[0] = x[0]
        out_x[1] = x[1] + (bbox[1] - bbox[4])
        out_x[2] = x[2]
        return out_x
   
    mpc = MultiPointConstraint(VF_space)
    mpc.create_periodic_constraint_geometrical(VF_space, locate_secondary, periodic_relation, bc_combined)
    mpc.finalize()

    mpc.add_constraint()

    print("MPC coeffs::", mpc.coefficients())
    # Define variational problem
    u_trial_function = ufl.TrialFunction(VF_space)
    v_test_function = ufl.TestFunction(VF_space)

    f = fem.Constant(msh, PETSc.ScalarType((0., 0., 0.)))

    Bil_form = fem.form(ufl.inner(sigma(u_trial_function), epsilon(v_test_function)) * ufl.dx(metadata={"quadrature_degree": ORDER}))
    Lin_form = fem.form(ufl.dot(f, v_test_function) * ufl.dx(metadata={"quadrature_degree": ORDER})) #+ ufl.dot(T, v) * ds

    # Setup MPC system

    problem = LinearProblem(Bil_form, Lin_form, mpc, bcs=bc_combined)

    solver = problem.solver

    # Give PETSc solver options a unique prefix
    solver_prefix = "dolfinx_mpc_solve_{}".format(id(solver))
    solver.setOptionsPrefix(solver_prefix)

    petsc_options: dict[str, Union[str, int, float]]
    petsc_options = {"ksp_type": "cg", "ksp_rtol": 1e-6, "pc_type": "hypre", "pc_hypre_type": "boomeramg",
                    "pc_hypre_boomeramg_max_iter": 1, "pc_hypre_boomeramg_cycle_type": "v"  # ,
                    # "pc_hypre_boomeramg_print_statistics": 1
                    }

    # Set PETSc options
    opts = PETSc.Options()  # type: ignore
    opts.prefixPush(solver_prefix)

    if petsc_options is not None:
        for k, v in petsc_options.items():
            opts[k] = v

    opts.prefixPop()
    solver.setFromOptions()

    uh = problem.solve()
    # solver.view()
    it = solver.getIterationNumber()
    print("Constrained solver iterations {0:d}".format(it))

    # Write solution to file
    outdir = Path("results_my_periodic")
    outdir.mkdir(exist_ok=True, parents=True)

    uh.name = "u_mpc"     

    with io.VTKFile(msh.comm, outdir / "demo_periodic_geometrical.vtk", "w") as vtk:
        vtk.write_mesh(msh)
        vtk.write_function(uh)
  

if __name__=="__main__":
    main();