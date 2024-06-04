using Gridap
using GridapGmsh
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.Visualization
using GridapPETSc
using GridapPETSc: PetscScalar, PetscInt
using SparseArrays
using SparseMatricesCSR
using IterativeSolvers


include("pbc_mpc_utils.jl")
include("penalty_mpc.jl")

input(filename::String) = "geometry/" * filename;
results(filename::String) = "results/" * filename;


tol = 1e-10
maxits = 1000
options = [
  "-ksp_type", "cg",
  "-ksp_monitor",
  "-ksp_rtol", "$tol",
  "-ksp_converged_reason",
  "-ksp_max_it", "$maxits",
  "-ksp_norm_type", "unpreconditioned",
  "-ksp_view",
  "-pc_type","gamg",
  "-pc_gamg_type","agg",
  "-mg_levels_esteig_ksp_type","cg",
  "-mg_coarse_sub_pc_type","cholesky",
  "-mg_coarse_sub_pc_factor_mat_ordering_type","nd",
  "-pc_gamg_process_eq_limit","50",
  "-pc_gamg_square_graph","0",
  "-pc_gamg_agg_nsmooths","1"];
const E = 70.0e9
const ν = 0.33
const λ = (E*ν)/((1+ν)*(1-2*ν))
const μ = E/(2*(1+ν))
# Safely open and close PETSc
GridapPETSc.with(args=options) do
   
    msh_file = "cell_057.msh" |> input
    # gmsh.initialize()
    model = GmshDiscreteModel(msh_file, renumber=false)

    vtk_file = "imported_mesh" |> input;
    writevtk(model, vtk_file);



    σ(ε) = λ*tr(ε)*one(ε) + 2*μ*ε

    order = 1
    degree = order*2
    Ω = Triangulation(model)


    dΩ = Measure(Ω,degree)

    reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},1)
    boundary_tags = [
      "n1_first", 
      "n1_second",
      # "n2_first", 
      # "n2_second",
      # "n3_first", 
      # "n3_second"
      ];
    boundary_masks = [
      (true,false,false), 
      (true,false,false), 
      # (false,true,false),
      # (false,true,false),
      # (false,false,true),
      # (false,false,true),
      ]
    test_space_V = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags = boundary_tags)

    g1(x) = VectorValue(x[1]*x[3],  0.0, 0.0)
    g2(x) = VectorValue(0.0,  0.0, 0.0)
    g3(x) = VectorValue(0.0,  0.0, -x[1]^2/2)

    boundaries = [
      g1, 
      g1, 
      # g2, 
      # g2, 
      # g3, 
      # g3,
    ]

    trial_space_U = TrialFESpace(test_space_V, boundaries) # Создаем пространство тестовых функций.

    Ω = Triangulation(model)
    dΩ = Measure(Ω,degree)


    bilinear_form_stiffness(u,v) = ∫( ε(v) ⊙ (σ∘ε(u)) )*dΩ
    linear_form_stiffness(v) = 0; #∫(dot(f, v))*dΩ

    solver = LinearFESolver(PETScLinearSolver())

    # Assembling on a Julia matrix
    # with the same data layout as petsc
    # (efficient use of PETScLinearSolver)
    Tm = SparseMatrixCSR{0,PetscScalar,PetscInt}
    Tv = Vector{PetscScalar}
    Tx = get_vector_type(U)
    assem = SparseMatrixAssembler(Tm, Tv, trial_space_U, test_space_V)
    op = AffineFEOperator(bilinear_form_stiffness, linear_form_stiffness, trial_space_U, test_space_V, assem)
    uh = solve(solver,op)

    # stiffness_matrix = get_matrix(operator_K)
    # rhs_vector = get_vector(operator_K)


    # du = cg(stiffness_matrix, rhs_vector)
    # rhs_vector
    ##
    res_file = "curvature_new" |> results
    # uh_lin = FEFunction(trial_space_U,du)
    writevtk(Ω, res_file ,  cellfields=["uh"=>uh, "stress"=>σ∘ε(uh)])
end