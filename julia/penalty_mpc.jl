# 1------2------3
# |      |      |
# |      |      |
# 4------5------6
# |      |      |
# |      |      |
# 7------8------9


abstract type MPC{T} end 

struct PenaltyMPC{T} <: MPC{T}
    penalty::T
    penalty_matrix::SparseMatrixCSC{T, Int64}
    penalty_forces::SparseVector{T, Int64}
end


function PenaltyMPC(pairs::Union{Base.Iterators.Zip, AbstractVector{<:Tuple{<:Integer, <:Integer}}}, 
    offsets::AbstractVector{<:Real}, 
    num_dofs::Integer, 
    penalty)
    @assert length(pairs) == size(offsets,1)
    B = spzeros(length(pairs), num_dofs);
    du = spzeros(length(pairs));

    for (row, (i,j)) in enumerate(pairs)
        B[row, i] = 1;
        B[row, j] = -1;
        if abs(offsets[row]) > 1e-8
            du[row] = offsets[row];
        end
    end 
    
    return PenaltyMPC(penalty, B' * B, B' * du)
end
mutable struct Dof2DofMap
    map::Vector{Tuple{Int64, Int64}}
    offsets::Vector{Float64}
    length::Int64
end

function Dof2DofMap(n::Int64)
    v_map = Vector{Tuple{Int64, Int64}}(undef, n);
    v_off = Vector{Float64}(undef, n);
    return Dof2DofMap(v_map, v_off, n);
end

function PenaltyMPC(dof_map::Dof2DofMap, num_dofs, penalty)
    return PenaltyMPC(dof_map.map, dof_map.offsets, num_dofs, penalty)
end


function applyMPC!(
    mK::AbstractMatrix{<:Real}, 
    vF::AbstractVector{<:Real},
    mpc::PenaltyMPC)
    mK += mpc.penalty_matrix * mpc.penalty;
    vF += mpc.penalty_forces * mpc.penalty;
end

