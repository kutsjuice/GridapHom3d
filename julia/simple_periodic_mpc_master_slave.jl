using WriteVTK
using LinearAlgebra
using SparseArrays
using LinearSolve
using IterativeSolvers
using PyPlot

remove(a, b) = a[findall(!in(b), a)]
findin(a, b) = a[findall(in(b), a)]


tcf(filename::String) = (@__DIR__) * "\\" * filename;

H = 1;
W = 1;
num = 51;

x = LinRange(0, W, num);
y = LinRange(0, H, num);

XX = reshape(repeat(x, inner = (1,num) )', num^2);
YY = reshape(repeat(y, inner = (1,num)), num^2);

points = [XX'; YY']

connections = Matrix{Int64}(undef, 3, (num-1)^2 * 2)
k = 1;
for i in 1:num-1
    for j in 1:num-1
        p1 = (i-1) * num + j;
        p2 = p1 + 1;
        p3 = i*num +j +1;
        p4 = p3-1;
        connections[:, k] .= [p1,p3,p2];
        connections[:, k+1] .= [p1,p4,p3];
        k+=2; 
    end
end

cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, connection) for connection in eachcol(connections)]

filename = "geo" |> tcf
vtk_grid(filename, points, cells) do vtk
end


Em = 5.0e10
Er = 2.1e11
ν_m = 0.45
ν_r = 0.3

# Model
flag = 1;

function buildelasticitymatrix(E, ν)
    k1 = ν/(1-ν);
    k2 = 0.5*(1-2ν)/(1-ν);
    k3 = (1-ν)/2;

    return E/(1-ν^2)* [
        1  ν  0;
        ν  1  0;
        0  0 k3;
    ]
end

mEm = buildelasticitymatrix(Em, ν_m);
mEr = buildelasticitymatrix(Er, ν_r);


r_bounds = [0.3 0.7]


function defineElasticityMatrix(pts)
    if (all(pts[1, :] .>= (r_bounds[1] - 1e-4)) 
        && all(pts[1, :] .<= (r_bounds[2] + 1e-4))
        && all(pts[2, :] .>= (r_bounds[1] - 1e-4))
        && all(pts[2, :] .<= (r_bounds[2] + 1e-4))
        )
        return mEr;
    else
        return mEm;
    end
end



function calcMKLocal(pts::Matrix{Float64})::Matrix{Float64}
    mB = Matrix{Float64}(undef, 3,6);
    ijk = [1,2,3];


    for _ in eachindex(pts)
        i,j,k = ijk;
        circshift!(ijk,1)
        # a = pts[1,j]*pts[2,k] - pts[1,k]*pts[2,j];
        b = pts[2,j] - pts[2,k];
        c = pts[1,k] - pts[1,j];
        mB[:, 2i-1:2i] .= [
            b 0; 
            0 c; 
            c b
            ];
    end
    Δ2 = det([
        1 pts[1,1] pts[2, 1];
        1 pts[1,2] pts[2, 2];
        1 pts[1,3] pts[2, 3];
        ]); 
    mB /= Δ2;
    mE = defineElasticityMatrix(pts)

    mKloc = mB' * mE * mB * Δ2/2;
    return mKloc;
end

function stress(pts::AbstractMatrix{<:Real}, disps::AbstractMatrix{<:Real})::Vector{Float64}
    mB = Matrix{Float64}(undef, 3,6);
    ijk = [1,2,3];
    for _ in eachindex(pts)
        i,j,k = ijk;
        circshift!(ijk,1)
        # a = pts[1,j]*pts[2,k] - pts[1,k]*pts[2,j];
        b = pts[2,j] - pts[2,k];
        c = pts[1,k] - pts[1,j];
        mB[:, 2i-1:2i] .= [
            b 0; 
            0 c; 
            c b
            ];
    end
    Δ2 = det([
        1 pts[1,1] pts[2, 1];
        1 pts[1,2] pts[2, 2];
        1 pts[1,3] pts[2, 3];
        ]); 
    mB /= Δ2;
    mE = defineElasticityMatrix(pts)
    σ = (mE * mB * reshape(disps, 6, 1))[:];
    return σ
end

if num < 20
    mK = zeros(num^2 * 2, num^2 * 2);
    vF = zeros(num^2 * 2);
    USING_SPARSE = false;
else 
    mK = spzeros(num^2 * 2, num^2 * 2);
    vF = spzeros(num^2 * 2);
    USING_SPARSE = true;
end

for el in eachcol(connections)
    mKloc = calcMKLocal(points[:, el]);
    inds = [2*el[1]-1, 2*el[1], 2*el[2]-1, 2*el[2], 2*el[3]-1, 2*el[3] ]
    mK[inds, inds] += mKloc;
end
rank(mK)


fixed = findall((points[1,:] .≈ 0.5) 
                .&& (points[2,:] .≈ 0.5))
# fixed = findall(points[1,:] .≈ 0)
# loaded = findall(points[1,:] .== 1)
loaded = Vector{Int64}([])

function applydirichelet!(
                mK::AbstractMatrix{<:Number}, 
                vF::AbstractVector{<:Number},
                indices::AbstractVector{<:Integer}, 
                dirsmask::BitVector, 
                func::Function, 
                points::AbstractMatrix{<:Number}
                )

    dims = size(dirsmask,1);
    dirs = findall(dirsmask);
    for ind in indices
        
        p = points[:, ind];
        dofsvalue = func(p);

        for dir in dirs
            dof = dims*(ind-1) + dir;            
            mK[dof, dof] = 1;
            vF[dof] = dofsvalue[dir];
            
            for k in eachindex(vF)
                if dof != k
                    if mK[dof, k] != 0
                        vF[k] -= dofsvalue[dir] * mK[dof, k];
                        mK[dof,k] = mK[k, dof] = 0;
                    end
                end
            end
        end 
    end

end

rank(mK)

f(p) = [0,0];
applydirichelet!(mK, vF, fixed, BitVector([1,1]), f, points)
mK
if USING_SPARSE
    droptol!(mK, 1e-3)
end
mK
rank(mK)

function findnewindices(old_subset, new_indices)
    new_subset = similar((old_subset));
    
    for i in eachindex(old_subset)
        new_subset[i] = findfirst(old_subset[i] .== new_indices)
    end
    return new_subset;
end

function applyMPC(  
            mK::AbstractMatrix{<:Number}, 
            vF::AbstractVector{<:Number}, 
            primary_dofs::AbstractVector{<:Integer}, 
            secondary_dofs::AbstractVector{<:Integer}, 
            koefs::AbstractVector{<:Number}, 
            offset::AbstractVector{<:Number}
        )

    @assert size(primary_dofs, 1) == size(secondary_dofs,1)
    @assert size(primary_dofs, 1) == size(koefs,  1)
    @assert size(primary_dofs, 1) == size(offset, 1)
    
    n = size(vF, 1);
    m = n - size(secondary_dofs, 1);
    uncommitted = remove(1:n, [primary_dofs;secondary_dofs]);
    uncommitted_and_primary = remove(1:n, secondary_dofs);

    new_primary = findnewindices(primary_dofs, uncommitted_and_primary)
    new_uncommitted = findall(in(uncommitted), uncommitted_and_primary);

    mT = spzeros(n,m);
    vG = spzeros(n);
    for (j,i) in enumerate(uncommitted_and_primary)
        mT[i,j] = 1;
    end
    for el in eachindex(primary_dofs)
        i = secondary_dofs[el];
        j = new_primary[el];
        mT[i,j] = koefs[el];
        vG[i] = offset[el];
    end
    vFnew = mT' * (vF - mK * vG) 
    mKnew = mT' * mK * mT;

    return (mKnew, vFnew, uncommitted_and_primary, new_primary, mT)
end

const EPS = 1e-6

isupper(p) = abs(p[2] - H) < EPS;
islower(p) = abs(p[2]) < EPS;
leftbound(p) = isapprox(p[1], 0; atol = 1e-3)
rightbound(p) = isapprox(p[1], 1; atol = 1e-3)



function finddofs(points::AbstractMatrix{<:Number}, condition::Function, number_of_directions::Integer)
    indices = findall(condition.(eachcol(points)));
    dofs = Vector{Int64}(undef, size(indices, 1) * number_of_directions);
    k = 0;
    for ind in indices
        for direction in 1:number_of_directions
            k += 1; 
            dofs[k] = number_of_directions*(ind-1) + direction;
        end
    end
    return dofs
end

upper_dofs = finddofs(points, isupper, 2); 
lower_dofs = finddofs(points, islower, 2);
left_dofs = finddofs(points, leftbound, 2);
right_dofs = finddofs(points, rightbound, 2);


upper_dofs' 
lower_dofs'
left_dofs'
right_dofs' 



left_dofs = remove(left_dofs, upper_dofs); left_dofs = remove(left_dofs, lower_dofs); left_dofs'
right_dofs = remove(right_dofs, lower_dofs); right_dofs = remove(right_dofs,  upper_dofs); right_dofs'
# upper_dofs'
# lower_dofs'


up_lo_offsets = zeros(size(lower_dofs))
for (i, dof) in enumerate(lower_dofs)
    if ((dof+1) % 2) == 0
        up_lo_offsets[i] = 0.1
    end
end

le_ri_offsets = zeros(size(left_dofs))

primary_dofs = [upper_dofs; left_dofs]
secondary_dofs = [lower_dofs; right_dofs]
offsets = [up_lo_offsets; le_ri_offsets]
up_lo_koeffs = ones(length(lower_dofs));
le_ri_koeffs = ones(length(left_dofs));
koeffs = ones(length(primary_dofs))

points[:,((primary_dofs.+1) .÷2)]
points[:,((secondary_dofs.+1) .÷2)]

offsets'


findin(primary_dofs, secondary_dofs)
secondary_dofs'
mKn,vFn, uncommitted_and_primary, new_primary, mT = applyMPC(mK, vF, primary_dofs, secondary_dofs, koeffs, offsets);
# mKn,vFn, uncommitted_and_primary, new_primary, mT = applyMPC(mK, vF, upper_dofs, lower_dofs, up_lo_koeffs, up_lo_offsets)

primary_dofs'
new_primary'


uncommitted_and_primary'


prob = LinearProblem(mKn, vFn);
sol = solve(prob)
δn = sol.u;


function restoredisps(
    δr::AbstractVector{<:Number},
    primary_dofs::AbstractVector{<:Integer},
    secondary_dofs::AbstractVector{<:Integer}, 
    koefs::AbstractVector{<:Number}, 
    offset::AbstractVector{<:Number}
    )::Vector{Float64}

    δ = Vector{Float64}(undef, length(δr) + length(secondary_dofs))
    uncommitted_and_primary = remove(1:length(δ), secondary_dofs); 


    δ[uncommitted_and_primary] = δr;
    δ[secondary_dofs] .= δ[primary_dofs] .* koefs + offset;
    return δ
end

δ = restoredisps(δn, primary_dofs, secondary_dofs, koeffs, offsets)
# δ = restoredisps(δn, upper_dofs, lower_dofs, up_lo_koeffs, up_lo_offsets)


# δ[uncommitted_and_primary] = δn;
# δ[lower_dofs] = δn[new_primary]
# rank(mK)
# δ = mK \ vF
stresses = Matrix{Float64}(undef, 3, length(cells))
for i in 1:size(connections,2)
    stresses[:, i] = stress(points[:, connections[:,i]], reshape(δ, 2, num^2)[:, connections[:,i]]);
end
# points
elasticity = Vector{Float64}(undef, length(cells));
for i in 1:size(connections,2)
    pts = points[:, connections[:,i]];
    if (all(pts[1, :] .>= (r_bounds[1] - 1e-4)) 
        && all(pts[1, :] .<= (r_bounds[2] + 1e-4))
        && all(pts[2, :] .>= (r_bounds[1] - 1e-4))
        && all(pts[2, :] .<= (r_bounds[2] + 1e-4)))
        elasticity[i] = Er;
           # println("!")
    else
        elasticity[i] = Em;
    end
end

vtk_grid(filename*"master-slave", points, cells) do vtk

    disps = zeros(3, num^2);
    disps[1:2,:] = reshape(δ, 2, num^2)
    vtk["disps"] = disps;
    vtk["stress"] = stresses;
    vtk["Young's modulus"] = elasticity
end