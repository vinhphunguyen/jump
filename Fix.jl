module Fix

using Printf
using LinearAlgebra
using Solid

abstract type FixBase end

struct EmptyFix <: FixBase
end

struct DisplacementFix <: FixBase
    solids::Vector{Solid2D}
    file
    index

    function DisplacementFix(solids,point,dir)
        index      = [-10 -10]
        println(point)
        for s = 1:length(solids)
            solid = solids[s]
            pos   = solid.pos
            for p = 1:solid.parCount
                x = pos[p][1]
                y = pos[p][2]
                if sqrt((x - point[1])^2 + (y - point[2])^2) < 0.9
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
    file
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
    file
    function StressFix(solids,filename)
        if (isfile(filename))
            rm(filename)
        end
        file = open(filename, "a")
        new(solids,file)
    end
end

function compute(fix::EmptyFix,time)
end

function closeFile(fix::EmptyFix)
end

function compute(fix::DisplacementFix,time)
    #println(fix.index)
    x = fix.solids[fix.index[1]].pos[fix.index[2]]
    write(fix.file, "$(time) $(x[1]) $(x[2])\n")
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

export FixBase, EmptyFix, EnergiesFix, DisplacementFix, StressFix, MaximumPlasticStrainFix
export compute, closeFile, compute_femp
end
