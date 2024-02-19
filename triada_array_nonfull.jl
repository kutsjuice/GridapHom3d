using Gmsh;
using Random;
using StaticArrays;
using StatsBase

struct TriadaArray
    bounds::SMatrix{3, 2, Float64}; #array with RVE bounds
    center::SVector{3, Float64}; #array with center of the RVE
    diam::Float64; #fibers' diameter
    intersection::Float64; #intersection ratio
    meshsize::Int64; # on each basic plane there would be rectangular grid with meshsize*meshsize fibers. So the overal number of fibers within the RVE will be 3*meshsize^2
    step::Float64; # distance between closest fibers in a row
    sfnum::Int64 # actual number of fibers on each basic plane. Range from 1 to meshsize^2
end

# function, which calculates cartesian indices (i,j) on rectangular grid by absolute index (k)
function abs2cartesian(inds::AbstractVector{<:Integer}, meshsize::Integer)::Vector{Tuple{Int64, Int64}}
    res_ = Vector{Tuple{Int64, Int64}}(undef, length(inds));
    for (k,ind) in enumerate(inds)
        i = ind ÷ meshsize + 1;
        j = ind % meshsize;
        if j == 0
            i -= 1;
            j = meshsize;
        end
        # a = (i,j);
        # display(typeof(a));
        res_[k] = (i,j);
    end
    return res_
end

#since we want to obtain a single body, all fibers should be intersected. to avoid floating fibers, i calculate the cartesian indices of fiber on third plane by indices on first and second plane
function thirdsideindices(inds1::Vector{Tuple{Int64, Int64}}, inds2::Vector{Tuple{Int64, Int64}}, inner::Bool = true)::Vector{Tuple{Int64, Int64}}
    res_ = Vector{Tuple{Int64, Int64}}(undef, length(inds1));

    for i in eachindex(inds1)
        if inner
            res_[i] = (inds1[i][1], inds2[i][2]);
        else
            res_[i] = (inds1[i][2], inds2[i][1]);
        end
    end
    return res_;
end


# main function, which build the geometry
function addTriadaArray(tr::TriadaArray, tag = -1)
    #ich fiber in the triad shifted from the center along both axes, perpendicular to fiber's axis
    offset = (1-tr.intersection)*tr.diam/2;

    #initialize fibers arrays of fibers tags
    I_s = Vector{Int64}(undef, tr.sfnum);
    J_s = Vector{Int64}(undef, tr.sfnum);
    K_s = Vector{Int64}(undef, tr.sfnum);

    # sample *sfnum* fibers from *meshsize^2*
    xy = sample(1:tr.meshsize^2, tr.sfnum; replace = false, ordered = true);
    yz = sample(1:tr.meshsize^2, tr.sfnum; replace = false, ordered = true);
    
    # calculate their cartesian indices
    XY_inds = abs2cartesian(xy, tr.meshsize);
    YZ_inds = abs2cartesian(yz, tr.meshsize);
    # caclulate indices for fibers on third plane
    XZ_inds = unique(thirdsideindices(XY_inds, YZ_inds));

    #it may happen, that for full bonding there would be enough lower number of fibers on third plane. In this case we need to add additional fibers in it. 
    while length(XZ_inds) < tr.sfnum
        rf = Tuple(rand(1:tr.meshsize, 2));
        if rf in XZ_inds
            continue;
        else
            push!(XZ_inds, rf)
        end
    end

    # generate fibers by cartesian indices in rectangular grid points on each plane
    for (k,(i, j)) in enumerate(XY_inds)
        x = tr.center[1] + offset + (i-1)*step - step*(tr.meshsize-1)/2;
        y = tr.center[2] + offset + (j-1)*step - step*(tr.meshsize-1)/2;

        zfiber = [  x  x;
                    y  y;
                    tr.bounds[3,1] tr.bounds[3,2];
                ];
        K_s[k] = gmsh.model.occ.addCylinder(zfiber[:, 1]..., (zfiber[:,2]-zfiber[:, 1])..., tr.diam/2);
    end

    for (k, (i, j)) in enumerate(YZ_inds)
        y = tr.center[2] - offset + (i-1)*step - step*(tr.meshsize-1)/2;
        z = tr.center[3] - offset + (j-1)*step - step*(tr.meshsize-1)/2;

        xfiber = [  tr.bounds[1,1] tr.bounds[1,2];
                    y  y;
                    z  z;
                ];
        
        I_s[k] = gmsh.model.occ.addCylinder(xfiber[:, 1]..., (xfiber[:,2]-xfiber[:, 1])..., tr.diam/2);
    
    end

    for (k, (i, j)) in enumerate(XZ_inds)

        x = tr.center[1] - offset + (i-1)*step - step*(tr.meshsize-1)/2;
        z = tr.center[3] + offset + (j-1)*step - step*(tr.meshsize-1)/2;

        yfiber = [  x x;
                    tr.bounds[2,1] tr.bounds[2,2];
                    z z;
                ];
    
        J_s[k] = gmsh.model.occ.addCylinder(yfiber[:, 1]..., (yfiber[:,2]-yfiber[:, 1])..., tr.diam/2);
    end


    # fuse all fibers to one body; we use I_s and J_s fibers as object
    objcDimTags = [(3, i) for i in vcat(I_s, J_s)];
    # and K_s as instrument for fusing
    toolDimTags = [(3, i) for i in K_s];

    # fusing
    bodies = gmsh.model.occ.fuse(objcDimTags, toolDimTags, tag);

    
    # return diamtags of fused body
    return bodies;
end

# setting initial parameters
diam = 25.5/1000
intersection = 0.25;


repnum = 6;
step = 50/1000 # 
# step = 0.549; # porosity 0.7
# step = 0.474; # porosity 0.6
# step = 0.423; # porosity 0.5
fibnum = 12;

width = step*repnum;
height = step*repnum;
depth = step*repnum;

bounds = SMatrix{3,2, Float64}([-width/2  width/2;
                                -height/2 height/2;
                                -depth/2  depth/2])

network = TriadaArray(bounds, [0,0,0], diam, intersection, repnum, step, fibnum)

# generation of the network. Since gmsh have to be finished anyway, I am using try...finally... structure. 
# I am recommending to use the similar in python for that, becouse when debugging, error may occure and we need gmsh to be finilized
try
    gmsh.initialize()

    # gmsh.option.setNumber("General.Terminal", 0);
    KK=12;
    # meshsize = 1000*(width^3)/(sqrt(2)^i)
    MSFC = 3;
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", MSFC);
    gmsh.option.setNumber("Mesh.MeshSizeMax", diam/2);
    gmsh.option.setNumber("Mesh.MeshSizeMin", diam/5);
    gmsh.option.setNumber("Mesh.MaxNumThreads3D", 8);
    gmsh.option.setNumber("Mesh.ElementOrder", 2);
    gmsh.option.setNumber("Mesh.Algorithm3D", 9);
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2);

    
    bodies = addTriadaArray(network);

    vₛ = gmsh.model.occ.getMass(bodies[1][1]...);
    vₜ = width*depth*height;
    
    por = 1 - vₛ/vₜ;
    # println("Porosity = $(1 - vₛ/vₜ)");

    # println("Bodies = ", bodies)        
    cutters = []
    push!(cutters, (3, gmsh.model.occ.addBox(4*bounds[1,1], 4*bounds[2,1], 4*bounds[3,1], 3*width/2, 3*height, 3*depth)));
    for i in 1:3
        push!(cutters, gmsh.model.occ.copy([cutters[end]])[1]);
        rot = gmsh.model.occ.rotate([cutters[end]], 0, 0, 0, 0, 0, 1, π/2);
    end
    push!(cutters, gmsh.model.occ.copy([cutters[end]])[1]);
    gmsh.model.occ.rotate(cutters[end], 0, 0, 0, 1, 0, 0, π/2);
    push!(cutters, gmsh.model.occ.copy([cutters[end]])[1]);
    gmsh.model.occ.rotate(cutters[end], 0, 0, 0, 1, 0, 0, π);

    gmsh.model.occ.cut(bodies[1][1], cutters)


    # gmsh.model.occ.cut(volumes, cutters)

    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)
    subdir="realizations/"
    filename = "POROSITY=$(round(por, digits=3))_DIAM=$(round(diam,digits=5))_FILLNES=$(round(fibnum/repnum^2, digits=3))_00";
    files = readdir(subdir)
    if filename*".msh" in files
        i=1;
        while filename[1:end-2]*lpad(string(i), 2, '0')*".msh" in files
            i+=1;
        end
        filename = filename[1:end-2]*lpad(string(i), 2, '0');
    end
    gmsh.write(subdir*"$(filename).msh")
    gmsh.write(subdir*"$(filename).vtk")

    # push!(outDimTags, triada_)
finally
    gmsh.finalize()
end

