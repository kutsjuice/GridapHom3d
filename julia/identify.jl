using Gridap

# prepare model and all related FE spaces
operator_K = AffineFEOperator(...)
operator_M = AffineFEOperator(...)
nodes = get_nodes_by_tag(model, tag);


springs = maximum(stiffness_matrix) * ones(length(nodes) * 3)


function calc_eigen_for_spring(springs)

    mass_matrix, stiffness_matrix = build_unconstrained_matrices(operator_K, operator_M);
    for (i,node) in enumerate(nodes)
        
        x_dof = node*3 -2
        stiffness = springs[3*i-2]
        y_dof = node*3 -1
        stiffness = springs[3*i-1]
        z_dof = node*3
        stiffness = springs[3*i]

        stiffness_matrix[x_dof, x_dof] += stiffness
        stiffness_matrix[y_dof, y_dof] += stiffness
        stiffness_matrix[z_dof, z_dof] += stiffness
    end
    v, frequences = calc_eigen_freq(stiffness_matrix, mass_matrix)
    return frequences;
end

function calc_jacobian(func, point)
    jac = 1 # compute jacobian of fun in point
    return jacobian
end

using ForwardDiff


A * x = f
x = inv(A) * f
x∈R^(N)
f∈R^(M)

A∈R^(M, N)

[A^T * A] * x = [A^T * f]
x = {[A^T * A]^(-1) * A^T} * f
A^+ = {[A^T * A]^(-1) * A^T}






# start iteration process from article
solved = false
for i in 1:1000
    error = calc_eigen_for_spring(s) - f0
    if norm(error) > EPS
        ds = calc_jacobian(calc_eigen_for_spring, spring) \ error;
        spring -= ds
    else
        solved = true
        break;
    end
end

if solved
    println("solution = $spring")
else
    println("reached 1000 itarations;")
end