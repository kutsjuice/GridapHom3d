using LinearAlgebra
using SparseArrays

# 1------2------3
# |      |      |
# |      |      |
# 4------5------6
# |      |      |
# |      |      |
# 7------8------9

B = zeros(6, 9)
d = zeros(6)
top2btm = [
    1 7;
    2 8;
    3 9;
]

lft2rht = [
    1 3;
    4 6;
    7 9;
]

function buildB(relations, num_dofs)
    B = zeros(size(relations, 1), num_dofs);
    for (row, (i,j)) in enumerate(eachrow(relations))
        B[row, i] = 1;
        B[row, j] = -1;
    end
    return B
end

B = buildB([top2btm;lft2rht], 9)

##

function makesparce!(L)
    for i in eachindex(L)
        if abs(L[i]) < 1e-6
            L[i] = 0.0
        end
    end
end
##
LT = nullspace(B)'

LT' * LT
# L = rref(L')'
# makesparce!(L)
# L
# L[:,3] = L[:,3] + L[:,1]

L = rref(LT)'
##

 L' * L

function rref(A::AbstractMatrix)
    n,m = size(A)
    ord_X = collect(1:m)
    
    A_pr = copy(A)          # делаем присоединённую матрицу
    
    for i in 1:n            # главный цикл - осуществляет операцию для каждой строки
        k = argmax(abs.(A_pr[i, i:end]))
        if maximum(abs.(A_pr[i, :])) < 1e-6
            continue
        end

        if i!=k
            A_pr[:, k+i-1], A_pr[:, i] = A_pr[:, i], A_pr[:, k+i-1]
            ord_X[i], ord_X[k+i-1] = ord_X[k+i-1], ord_X[i]
        end
        A_pr[:, i] = A_pr[:, i] ./ A_pr[i,i]
        for j in i+1:n      # вложенный цикл для вычитания
            for k in i:m
                A_pr[j, k] -= A_pr[i, k] * A_pr[j, i] # из j-й строки вычесть i-ю, умноженную на соответствующие коэффициенты. Чтобы сделать это со всеми столбацами, используется срез.
            end
        end
    end

    # обратны ход
    for i in 2:n
        for j in 1:i-1
            for k in j:m
                A_pr[j,k] -= A_pr[i,k] * A_pr[j,i] 
            end
        end
    end

    return A_pr[:, ord_X]
end

ge(L)

B = nullspace(L)
BT = ge(B')

BT[1, :] /= BT[1, 1]
BT[2, :] /= BT[2, 2]
BT[3, :] /= BT[3, 3]
BT'

# 1------2------3
# |      |      |
# |      |      |
# 4------5------6
# |      |      |
# |      |      |
# 7------8------9
# u1 = 1 * u7 + 0; - (1, 7, 1.0, 0.0)
# u2 = u8; 
# u3 = u9 
# u1 = u3; 
# u4 = u6; 
# u7 = u9
# 
# 
#

##
top2btm = [
    1 7;
    2 8;
    3 9;
]

lft2rht = [
    1 3;
    4 6;
    7 9;
]


# mutable struct SingleMPC
#     master::Int64
#     slave::Int64
#     coefficient::Float64
#     offset::Float64
# end

mutable struct SingleMPCSet
    masters::Dict{Int64, Vector{Int64}}
    slaves::Dict{Int64, Int64}
    mpc::Vector{Tuple{Int64, Int64, Float64, Float64}}
end


# SingleMPC(master::Int64, slave::Int64; coefficient::Float64 = 1.0, offset::Float64 = 0.0) = SingleMPC(master, slave, coefficient, offset);

#             | master | slave |
#             |        |       |
# set.masters | 1 | 0  | 1 | 0 |
#             |        |       |
# set.slaves  | 1 | 0  | 1 | 0 |
#
# in overal we have 9 cases - [m ∈ [M, S, ∅]] × [s ∈ [M, S, ∅]]
# 1. m∈M & s∈M - remap M[s] mpc on new m
# 2. m∈M & s∈S - check, if the new mpc is the same
# 3. m∈M & s∈∅ - just add new mpc
# 4. m∈S & s∈M - invert mpc and go to p.2.
# 5. m∈S & s∈S - remap m on um from S[m] and go to p.2 
# 6. m∈S & s∈∅ - remap m on um from S[m] and go to p.3
# 7. m∈∅ & s∈M - remap M[s] mpc on new m 
# 8. m∈∅ & s∈S - go to p.2
# 9. m∈∅ & s∈∅ - go to p.3
#
if master ∈ keys(set.masters)
    if slave ∈ keys(set.slaves)
        
    elseif slave ∈ keys(set.masters)
        
    else
        
    end
elseif master ∈ keys(set.slaves)
    if slave ∈ keys(set.slaves)

    elseif slave ∈ keys(set.masters)

    else

    end 

else
    if slave ∈ keys(set.slaves)

    elseif slave ∈ keys(set.masters)

    else

    end
end

s = m*c + o 
# function addmpc!(set::SingleMPCSet, master, slave; coeff, offset)
#     if haskey(set.slaves, master)
#         if !haskey(set.masters, slave)
#             u_mpc = set.mpc[set.slaves[master]];
#             u_master = u_mpc[1];
#             u_coeff = u_mpc[3];
#             u_offset = u_mpc[4];
#             # usl = ucf*ums + uof
#             # sl = cf*ms + of = cf*(ucf*ums + uof) + of = (cf*ucf)*ums + (cf*uof + of)
#             push!(set.mpc, (u_master, slave, u_coeff*coeff, coeff*u_offset + offset))
#             push!(set.masters[u_master, length(set.mpc)])
#             set.slaves[slave] = length(set.mpc);
#         else

#         end
#     elseif 

#     end
# end

for (ms, sl) in top2btm
    addmpc!()
end







##

