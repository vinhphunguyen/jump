# ----------------------------------------------------------------------
#
#                    ***       JUMP       ***
#                Material Point Method in Julia
#
# Copyright (2020) Vinh Phu Nguyen, phu.nguyen@monash.edu
# Civil Engineering, Monash University
# Clayton VIC 3800, Australia
# This software is distributed under the GNU General Public License.
#
# -----------------------------------------------------------------------

module Output
import Glob
import PyPlot

using Solid
using Material

abstract type OutputType end

struct PyPlotOutput <: OutputType
	interval      ::Int64
	dir           ::String
	figTitle      ::String
    figSize       ::Tuple{Float64,Float64}

	function PyPlotOutput(interval::Int64,dir::String,figTitle::String,figSize::Tuple{Float64,Float64})
		if !isdir(dir)	mkdir(dir) end
        new(interval,dir,figTitle,figSize)
    end
end

struct OvitoOutput <: OutputType
	interval      ::Int64
	dir           ::String
	outs          ::Vector{String}

	function OvitoOutput(interval::Int64,dir::String,outs)
		if isdir(dir)
			#command = 'rm *.LAMMPS'
			dumpfiles = Glob.glob(string(dir,"*.LAMMPS"))
			if (length(dumpfiles) > 0 )
				[rm(file)  for file in dumpfiles]
			end
		else
			mkdir(dir)
		end
        new(interval,dir,outs)
    end
end

function plotParticles(plot::PyPlotOutput,solids::Vector{Solid2D},
                       lims::Vector{Float64},ncells::Vector{Int64},counter::Int64)
    fileName = string(plot.dir,"$(Int(counter)).pdf")
	#array_x         = [toXArray(solids[i]) for i = 1:length(solids)]
	#array_y         = [toYArray(solids[i]) for i = 1:length(solids)]
	array_x         = Vector{Float64}(undef,0)
	array_y         = Vector{Float64}(undef,0)
	@inbounds for s=1:length(solids)
		solid=solids[s]
		for p=1:solid.parCount
			push!(array_x,solid.pos[p][1])
			push!(array_y,solid.pos[p][2])
		end
	end

	#    array_color     = Array{Real}(iMaterialPoints, 3)
	# array_size      = Array{Real}(iMaterialPoints, 1)
	#    for iIndex in 1:1:iMaterialPoints
	#      array_color[iIndex, :] = [1.0, 0.0, 0.0]#[thisGrid.GridPoints[iIndex].fMass/iMaterialPoints, 0.0, 0.0]
	#   array_size[iIndex, :] = [5.0]
	#    end

    xlim = lims[1]
    ylim = lims[2]
    dxx  = xlim/ncells[1]
    dyy  = ylim/ncells[2]
	pyPlot01 = PyPlot.gca()
	# pyPlot01 = PyPlot.subplot2grid((1,1), (0,0), colspan=1, rowspan=1, aspect="equal")
	PyPlot.scatter(array_x, array_y, s=1)
	pyPlot01[:spines]["top"][:set_color]("gray")
	pyPlot01[:spines]["right"][:set_color]("gray")
	pyPlot01[:spines]["bottom"][:set_color]("gray")
	pyPlot01[:spines]["left"][:set_color]("gray")
	# pyPlot01[:axhline](linewidth=4, color="g")
	# pyPlot01[:axvline](linewidth=4, color="g")
	pyPlot01[:set_xlim](0.0, xlim)
	pyPlot01[:set_ylim](0.0, ylim)
	# pyPlot01[:set_xlabel]("")
	# pyPlot01[:set_ylabel]("")
	pyPlot01[:grid](b=true, which="both", color="gray", linestyle="-", linewidth=0.5)
	pyPlot01[:set_axisbelow](true)
	pyPlot01[:set_xticks]([])# empty to have no major ticks and grids
	pyPlot01[:set_xticks](collect(0.0:dxx:xlim),minor=true)
	pyPlot01[:set_yticks]([])# empty to have no major ticks and grids
	pyPlot01[:set_yticks](collect(0.0:dyy:ylim),minor=true)

	#PyPlot.show()
	PyPlot.pause(0.05)
	PyPlot.savefig(fileName, bbox_inches="tight")
	PyPlot.clf()
	PyPlot.gcf()
end

#ITEM: TIMESTEP
# 0
# ITEM: NUMBER OF ATOMS
# 334328
# ITEM: BOX BOUNDS sm sm sm
# -17 50
# -25 25
# -25 25
#ITEM: ATOMS id type x y z damage s11 s22 s33 s12 s13 s23 vx T ienergy
#0 1 10.5401 -1.50408 -0.647102 0 -16.9889 -19.2388 -13.7454 11.4435 4.68828 -2.59562 3346.25 298 0
function plotParticles(plot::OvitoOutput,solids,
          lims::Vector{Float64},ncells::Vector{Int64},counter::Int64)
	parCount = 0
	for s=1:length(solids)
		parCount += solids[s].parCount
	end
	fileName = string(plot.dir,"dump_p.","$(Int(counter)).LAMMPS")
	println(fileName)
	file = open(fileName, "a")
    write(file, "ITEM: TIMESTEP \n")
    write(file, "0\n")
    write(file, "ITEM: NUMBER OF ATOMS\n")
	write(file, "$(parCount) \n")
	write(file, "ITEM: BOX BOUNDS sm sm sm\n")
	write(file, "0 $(lims[1])\n")
	write(file, "0 $(lims[2])\n")
	write(file, "0 0\n")		 # 3D not work with this

	write(file, "ITEM: ATOMS id type x y z ")
	@inbounds for v in plot.outs
		write(file, v)
		write(file, " ")
	end
	write(file, "\n")
    id = 0
	@inbounds for s=1:length(solids)
		x       = solids[s].pos
		elastic = false
		mat     = solids[s].mat
		stress  = solids[s].stress
		velocity = solids[s].velocity
		if ( typeof(mat) <: ElasticMaterial ) || ( typeof(mat) <: RigidMaterial )
			elastic = true
		end
		@inbounds for p = 1:length(x)
			id += 1
			write(file, "$(id) $(s) $(x[p][1]) $(x[p][2]) 0 ")
			for v in plot.outs
				if v == "pstrain"
					if elastic
						write(file, "0 ")
					else
				        write(file, "$(mat.alpha[p]) ")
					end
				elseif v == "vonMises"
					if elastic
						write(file, "0 ")
					else
				        write(file, "$(mat.vmStr[p]) ")
					end
				elseif v == "sigmaxx"
					write(file, "$(stress[p][1,1]) ")
				elseif v == "sigmayy"
					write(file, "$(stress[p][2,2]) ")
				elseif v == "pressure"
					write(file, "$(stress[p][1,1]+
					               stress[p][2,2]) ")
				elseif v == "vx"
					write(file, "$(velocity[p][1])  ")
				elseif v == "vy"
					write(file, "$(velocity[p][2])  ")
				elseif v == "damage"
					write(file, "$(solids[s].damage[p])  ")
				elseif v == "crackforce"
					write(file, "$(solids[s].crackForce[p])  ")
                elseif v == "color"
				    write(file, "$(solids[s].color[p])  ")
				end
			end
			write(file, "\n")
		end
	end
	close(file)
end

function plotParticles(plot::OvitoOutput,solids::Vector{Solid3D},
          lims::Vector{Float64},ncells::Vector{Int64},counter::Int64)
	parCount = 0
	for s=1:length(solids)
		parCount += solids[s].parCount
	end
	fileName = string(plot.dir,"dump_p.","$(Int(counter)).LAMMPS")
	println(fileName)
	file = open(fileName, "a")
    write(file, "ITEM: TIMESTEP \n")
    write(file, "0\n")
    write(file, "ITEM: NUMBER OF ATOMS\n")
	write(file, "$(parCount) \n")
	write(file, "ITEM: BOX BOUNDS sm sm sm\n")
	write(file, "0 $(lims[1])\n")
	write(file, "0 $(lims[2])\n")
	write(file, "0 $(lims[3])\n")

	write(file, "ITEM: ATOMS id type x y z")
	@inbounds for v in plot.outs
		write(file, v)
		write(file, " ")
	end
	write(file, "\n")
    id = 0
	@inbounds for s=1:length(solids)
		x       = solids[s].pos
		elastic = false
		mat     = solids[s].mat
		if ( typeof(mat) <: ElasticMaterial ) elastic = true end
		@inbounds for p = 1:length(x)
			id += 1
			write(file, "$(id) $(s) $(x[p][1]) $(x[p][2]) $(x[p][3] ) ")
			for v in plot.outs
				if v == "pstrain"
					if elastic
						write(file, "0 ")
					else
				        write(file, "$(mat.alpha[p]) ")
					end
				elseif v == "vonMises"
					if elastic
						write(file, "0 ")
					else
				        write(file, "$(mat.vmStr[p]) ")
					end
				elseif v == "vx"
					write(file, "$(solids[s].velocity[p][1])  ")
				elseif v == "vy"
					write(file, "$(solids[s].velocity[p][2])  ")
				elseif v == "vz"
					write(file, "$(solids[s].velocity[p][3])  ")
				# elseif v == "surf"
				# 	println(solids[s].color[p])
				# 	write(file, "$(solids[s].color[p])  ")
				elseif v == "pressure"
					write(file, "$(solids[s].stress[p][1,1]+
							   	   solids[s].stress[p][2,2]+solids[s].stress[p][3,3])  ")
				end
			end
			write(file, "\n")
		end
	end
	close(file)
end

function plotParticles(plot::OvitoOutput,grid,counter::Int64)
	fileName = string(plot.dir,"dump_g.","$(Int(counter)).LAMMPS")
	println(fileName)
	file = open(fileName, "a")
    write(file, "ITEM: TIMESTEP \n")
    write(file, "0\n")
    write(file, "ITEM: NUMBER OF ATOMS\n")
	write(file, "$(grid.nodeCount) \n")
	write(file, "ITEM: BOX BOUNDS sm sm sm\n")
	write(file, "0 $(grid.lx)\n")
	write(file, "0 $(grid.ly)\n")
	write(file, "0 0\n")		 # 3D not work with this

	write(file, "ITEM: ATOMS id type x y z vx vy \n")
    xx=grid.pos
	vel = grid.momentum
	@inbounds for i = 1:grid.nodeCount
		write(file, "$(i) 1 $(xx[i][1]) $(xx[i][2]) 0 $(vel[i][1]) $(vel[i][2]) \n")
	end
	close(file)
end

export OutputType, PyPlotOutput, OvitoOutput
export plotParticles, writeParticles

end
