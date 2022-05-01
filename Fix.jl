module Fix

using Printf
using LinearAlgebra
using StaticArrays

using Solid
using Fem
using Grid
using Mesh


abstract type FixBase end

struct EmptyFix <: FixBase
end

struct DisplacementFix <: FixBase
    solids#::Vector{Solid3D}
    file  ::IOStream
    index ::Matrix{Int64}

    function DisplacementFix(solids,point::SVector{2,Float64},dir)
        index      = [-10 -10]
        println(point)
        for s = 1:length(solids)
            solid = solids[s]
            pos   = solid.pos
            for p = 1:solid.parCount
                x = pos[p][1]
                y = pos[p][2]
                if sqrt((x - point[1])^2 + (y - point[2])^2) < 1e-10
                    index = [s p]
                    #println(index)
                    break
                end
            end
        end

        if (index[2] < 0 ) error("No particle with that position found\n") end

        filename = string(dir,"recorded-position.txt")

        if (isfile(filename))
            rm(filename)
        end

        file = open(filename, "a")

        new(solids,file,index)
    end

    function DisplacementFix(solids,point::SVector{3,Float64},dir)
        index      = [-10 -10]
        println(point)
        for s = 1:length(solids)
            solid = solids[s]
            pos   = solid.pos
            for p = 1:solid.parCount
                x = pos[p][1]
                y = pos[p][2]
                z = pos[p][3]
                if sqrt((x - point[1])^2 + (y - point[2])^2 + (z - point[3])^2) < 1e-10
                    index = [s p]
                    #println(index)
                    break
                end
            end
        end

        if (index[2] < 0 ) error("No particle with that position found\n") end

        filename = string(dir,"recorded-position.txt")

        if (isfile(filename))
            rm(filename)
        end

        file = open(filename, "a")

        new(solids,file,index)
    end
end

struct EnergiesFix <: FixBase
    solids
    file::IOStream
	kinEnergy ::Vector{Float64}
	strEnergy ::Vector{Float64}
	recordTime::Vector{Float64}
	function EnergiesFix(solids,filename)
        if (isfile(filename))
            rm(filename)
        end
        file = open(filename, "a")
        new(solids,file,Vector{Float64}(),Vector{Float64}(),Vector{Float64}())
    end
end

struct StressFix <: FixBase
    solids
    file::IOStream
    function StressFix(solids,filename)
        if (isfile(filename))
            rm(filename)
        end
        file = open(filename, "a")
        new(solids,file)
    end
end

struct CenterOfMassFix <: FixBase
    solid
    file::IOStream
    function CenterOfMassFix(solid,filename)
        if (isfile(filename))
            rm(filename)
        end
        file = open(filename, "a")
        new(solid,file)
    end
end

struct DisplacementFemFix <: FixBase
    solid
    file::IOStream
    nodeID::Int64

    function DisplacementFemFix(solid,dir,nodeID)
        filename = string(dir,"recorded-position.txt")
        if (isfile(filename))
            rm(filename)
        end

        file = open(filename, "a")

        new(solid,file,nodeID)
    end
end

struct ReactionForceFix <: FixBase
    solid::FEM3D
    file::IOStream
    function ReactionForceFix(solid,tag,filename)
        if (isfile(filename))
            rm(filename)
        end
        file = open(filename, "a")
        if haskey(solid.mesh.node_sets,tag ) == false 
           Mesh.create_node_set_from_element_set!(solid.mesh, tag)
        end
        new(solid,file)
    end
end

struct SpatialReactionForceFix <: FixBase
    solid::FEM3D
    file::IOStream
    tag::String
    function SpatialReactionForceFix(solid,tag,filename)
        if (isfile(filename))
            rm(filename)
        end
        file = open(filename, "a")
        if haskey(solid.mesh.node_sets,tag ) == false 
           Mesh.create_node_set_from_element_set!(solid.mesh, tag)
        end
        new(solid,file,tag)
    end
end


struct ReactionForceGridFix <: FixBase
    grid::Grid3D
    file::IOStream
    function ReactionForceGridFix(grid,filename)
        if (isfile(filename))
            rm(filename)
        end
        file = open(filename, "a")
        new(grid,file)
    end
end

###########################################
# define functions compute(...) here
###########################################

function compute(fix::EmptyFix,time)
end

function closeFile(fix::EmptyFix)
end

function compute(fix::DisplacementFix,time)
    #println(fix.index)
    x = fix.solids[fix.index[1]].pos[fix.index[2]]
    write(fix.file, "$(time) $(x[1]) $(x[2])\n")
end

function compute_femp(fix::EmptyFix,time)
end

function compute_femp(fix::DisplacementFemFix,time)
    #println(fix.index)
    x = fix.solid.pos[fix.nodeID][1] - fix.solid.pos0[fix.nodeID][1]
    y = fix.solid.pos[fix.nodeID][2] - fix.solid.pos0[fix.nodeID][2]
    z = 0
    if typeof(fix.solid) <: FEM3D z = fix.solid.pos[fix.nodeID][3] - fix.solid.pos0[fix.nodeID][3] end
    write(fix.file, "$(time) $(x) $(y) $(z)\n")
end

function compute_femp(fix::CenterOfMassFix,time)
    tmass = 0.
    xcm   = @SVector[0,0]
    mm    = fix.solid.mass
    for ip = 1:fix.solid.nodeCount
        tmass += mm[ip] 
        xcm   += mm[ip] * fix.solid.pos[ip]    
    end
    #if typeof(fix.solid) <: FEM3D z = fix.solid.pos[fix.nodeID][3] - fix.solid.pos0[fix.nodeID][3] end
    write(fix.file, "$(time) $(xcm[1]/tmass) $(xcm[2]/tmass)\n")
end

function compute(fix::EnergiesFix,time)
	KE     = 0.
	SE     = 0.
	solids = fix.solids
	for s = 1:length(solids)
		solid  = solids[s]

		mm     = solid.mass
		vv     = solid.velocity
		vol    = solid.volume
		stress = solid.stress
		strain = solid.strain

        for ip = 1:solid.parCount
          KE += 0.5 * mm[ip]  * dot(vv[ip],vv[ip])
          SE += 0.5 * vol[ip] * sum(stress[ip].*strain[ip])
        end
	end
	push!(fix.kinEnergy,KE)
	push!(fix.strEnergy,SE)
	push!(fix.recordTime,time)
    write(fix.file, "$(time) $(KE) $(SE)\n")
end

function compute_femp(fix::EnergiesFix,time)
    KE     = 0.
    SE     = 0.
    solids = fix.solids
    for s = 1:length(solids)
        solid  = solids[s]

        mm     = solid.mass
        vv     = solid.velocity
        vol    = solid.volume
        stress = solid.stress
        strain = solid.strain
        

        for ip = 1:solid.nodeCount
          KE += 0.5 * mm[ip]  * dot(vv[ip],vv[ip])
          #SE += 0.5 * vol[ip] * sum(stress[ip].*strain[ip])
        end

        for ip = 1:solid.parCount
          #KE += 0.5 * mm[ip]  * dot(vv[ip],vv[ip])
          SE += 0.5 * vol[ip] * sum(stress[ip].*strain[ip])
        end
    end
    push!(fix.kinEnergy,KE)
    push!(fix.strEnergy,SE)
    push!(fix.recordTime,time)
    write(fix.file, "$(time) $(KE) $(SE)\n")
end

function compute_femp(fix::ReactionForceFix,time)
   
    ids  = collect(fix.solid.mesh.node_sets["bottom"])
    force = fix.solid.fint

    fx = fy = fz = 0.
    for id in ids
        fx += force[id][1]
        fy += force[id][2]
        fz += force[id][3]
    end
    write(fix.file, "$(time) $(fx) $(fy) $(fz)\n")


    #write(fix.file, "$(time) $(fix.solid.reaction_forces[1]) $(fix.solid.reaction_forces[2]) $(fix.solid.reaction_forces[3])\n")
end


function compute_femp(fix::ReactionForceGridFix,time)
    grid = fix.grid
    force = grid.force
    fx = fy = fz = 0.
    for id in grid.bottomNodes
        fx += force[id][1]
        fy += force[id][2]
        fz += force[id][3]
    end
    write(fix.file, "$(time) $(fx) $(fy) $(fz)\n")
end

function compute_femp(fix::SpatialReactionForceFix,time)
    ids   = collect(fix.solid.mesh.node_sets[fix.tag])
    force = fix.solid.fint
    xx    = fix.solid.pos

    for id in ids
        x  = xx[id][1]
        fx = force[id][1]
        fy = force[id][2]
        fz = force[id][3]
        write(fix.file, "$(x) $(fx) $(fy) $(fz)\n")
    end        
end

function compute(fix::StressFix,time)
    #println(fix.index)    
    write(fix.file, "$(fix.solids[1].strain[1][1,1]) $(fix.solids[1].stress[1][1,1]) \n")
end

function closeFile(fix::FixBase)
    close(fix.file)
end

# function closeFile(fix::EnergiesFix)
#     close(fix.file)
# end

struct MaximumPlasticStrainFix <: FixBase

end

export FixBase, EmptyFix, EnergiesFix, DisplacementFix, StressFix, MaximumPlasticStrainFix, DisplacementFemFix, CenterOfMassFix
export ReactionForceFix, ReactionForceGridFix, SpatialReactionForceFix
export compute, closeFile, compute_femp
end
