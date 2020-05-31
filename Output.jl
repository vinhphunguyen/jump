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
using  WriteVTK
using  StaticArrays

using Solid
using Material
using Grid

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

struct VTKOutput <: OutputType
	interval      ::Int64
	dir           ::String
	outs          ::Vector{String}


	function VTKOutput(interval::Int64,dir::String,outs)
		if isdir(dir)
			#command = 'rm *.LAMMPS'
			dumpfiles = Glob.glob(string(dir,"*.vtu"))
			if (length(dumpfiles) > 0 )
				[rm(file)  for file in dumpfiles]
			end
		else
			mkdir(dir)
		end
        new(interval,dir,outs)
    end
end

function plotParticles(plot::PyPlotOutput,solids,
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
function plotParticles_2D(plot::OvitoOutput,solids,
          lims::Vector{Float64},ncells::Vector{Int64},counter::Int64)
	parCount = 0
	for s=1:length(solids)
		parCount += solids[s].parCount
	end
	fileName = string(plot.dir,"dump_p.","$(Int(counter)).LAMMPS")
	#println(fileName)
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

function plotParticles_3D(plot::OvitoOutput,solids,
          lims::Vector{Float64},ncells::Vector{Int64},counter::Int64) where {T<:MaterialType}
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

function plotGrid(plot::OvitoOutput,grid::Grid2D,counter::Int64)
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

function plotParticles_2D(plot::VTKOutput,solids,counter::Int64)
	my_vtk_file = string(plot.dir,"particle_","$(Int(counter))")
	nodeCount = 0
	for s=1:length(solids)
		nodeCount += solids[s].nodeCount
	end
	points      = zeros(2,nodeCount)
	cc = 1
	cells = MeshCell[]
	p     = Vector{Float64}(undef,0)
	velo  = zeros(2,nodeCount)
	for s=1:length(solids)
		solid = solids[s]
		xx    = solid.pos
		elems = solid.elems
		stress= solid.stress
		ve    = solid.velocity
		shift = (s-1)*solid.nodeCount
		for ip=1:solid.nodeCount
			points[1,cc] = xx[ip][1]
			points[2,cc] = xx[ip][2]
			velo[1,cc]   = ve[ip][1]
			velo[2,cc]   = ve[ip][2]
			cc += 1
		end
		for e=1:solid.parCount
			inds =elems[e,:] .+ shift
		    c    = MeshCell(VTKCellTypes.VTK_QUAD, inds)
            push!(cells, c)
            sigma = stress[e]
            push!(p, sigma[1,1]+sigma[2,2])
        end
	end
	vtkfile     = vtk_grid(my_vtk_file, points, cells)
	# write data 
	vtkfile["Pressure", VTKCellData()]  = p
	vtkfile["Velocity", VTKPointData()] = velo
	outfiles    = vtk_save(vtkfile)
end

function plotGrid(plot::VTKOutput,grid::Grid2D)
	my_vtk_file = string(plot.dir,"grid")

	points      = zeros(2,grid.nodeCount)

	xx=grid.pos
	
	@inbounds for i = 1:grid.nodeCount
		points[1,i] = xx[i][1] 
		points[2,i] = xx[i][2]
	end

	
	cells = MeshCell[]
	for j = 1:grid.nodeCountY-1
		for i = 1:grid.nodeCountX-1
			node1   = i + (j-1)*grid.nodeCountX
			node2   = node1 + 1
			node3   = node1 + grid.nodeCountX
			node4   = node3 + 1
            inds    = @SVector[node1,node2,node4,node3]
		    c    = MeshCell(VTKCellTypes.VTK_QUAD, inds)
            push!(cells, c)
        end
	end
	vtkfile     = vtk_grid(my_vtk_file, points, cells)
	outfiles    = vtk_save(vtkfile)
end

function plotGrid(plot::VTKOutput,grid::Grid3D)
	my_vtk_file = string(plot.dir,"grid")

	points      = zeros(3,grid.nodeCount)

	xx=grid.pos
	
	@inbounds for i = 1:grid.nodeCount
		points[1,i] = xx[i][1] 
		points[2,i] = xx[i][2]
		points[3,i] = xx[i][3]
	end

	
	cells = MeshCell[]
	noXY  = nodeCountX * nodeCountY
	for k = 1:grid.nodeCountZ-1
		for j = 1:grid.nodeCountY-1
			for i = 1:grid.nodeCountX-1
				node1   = i + (j-1)*grid.nodeCountX + (k-1)*noXY
				node2   = node1 + 1
				node3   = node1 + grid.nodeCountX
				node4   = node3 + 1

				node5   = node1 + noXY
				node6   = node2 + noXY
				node7   = node3 + noXY
				node8   = node4 + noXY

	            inds    = @SVector[node1,node2,node4,node3,node5,node6,node8,node7]
			    c       = MeshCell(VTKCellTypes.VTK_HEX, inds)
	            push!(cells, c)
	        end
		end
    end
	vtkfile     = vtk_grid(my_vtk_file, points, cells)
	outfiles    = vtk_save(vtkfile)
end

export OutputType, PyPlotOutput, OvitoOutput, VTKOutput
export plotGrid,plotParticles_2D, plotParticles_3D, writeParticles, plotGrid

end
