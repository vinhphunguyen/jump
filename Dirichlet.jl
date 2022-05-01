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

module Dirichlet

# This module implements functions for applying Dirichlet boundary conditions on a Eulerian grid
# data["dirichlet_grid"] = [("bottom",(0,1,0)),      # => fix bottom nodes on Y dir
#                           ("top",(1,0,0)),      # => fix top nodes on X dir
#                          ]

using StaticArrays
using Mesh
using Util
using Grid

# input file: data["dirichlet_grid"] = [("bottom",(0,1,0))] # => fix bottom nodes on Y dir
# 3D grid 

function fix_Dirichlet_grid(grid::Grid3D,data;ghostcell=false)
 	if haskey(data, "dirichlet_grid") == false return end 

    dirichlet_grid = data["dirichlet_grid"]

 	for (bnd,fix) in dirichlet_grid

 		if     bnd == "bottom" grid.fixedNodes[:,grid.bottomNodes] .= fix
 		elseif bnd == "top"    grid.fixedNodes[:,grid.topNodes]    .= fix
 		elseif bnd == "left"   grid.fixedNodes[:,grid.leftNodes]   .= fix
 		elseif bnd == "right"  grid.fixedNodes[:,grid.rightNodes]  .= fix
 		elseif bnd == "front"  grid.fixedNodes[:,grid.frontNodes]  .= fix
 		elseif bnd == "back"   grid.fixedNodes[:,grid.backNodes]   .= fix
 		else  error("Not recognized grid faces!!!")
 		end
 	end
end

# input file: data["dirichlet_grid"] = [("bottom",(0,1,0))] # => fix bottom nodes on Y dir
# NOTE: ghostcell=true, only left edge is done!!!

function fix_Dirichlet_grid(grid::Grid2D,data;ghostcell=false)
 	if haskey(data, "dirichlet_grid") == false return end 

    dirichlet_grid = data["dirichlet_grid"]

    # grid.fixedNodes = a matrix of size (2,nodeCount)
 	for (bnd,fix) in dirichlet_grid

 		if     bnd == "bottom" 
 		   grid.fixedNodes[:,grid.bottomNodes] .= fix

           if ghostcell == true
 		   	  grid.fixedNodes[:,grid.bottomNodes .+ grid.nodeCountX] .= fix	
 		   end 
 		elseif bnd == "top"    
 			grid.fixedNodes[:,grid.topNodes]    .= fix

 			if ghostcell == true
 			   grid.fixedNodes[:,grid.topNodes .- grid.nodeCountX]    .= fix
 			end
 		elseif bnd == "left"   
            grid.fixedNodes[:,grid.leftNodes]   .= fix
            if ghostcell == true
           	   grid.fixedNodes[:,grid.leftNodes.+1]   .= fix
            end
 		elseif bnd == "right"  
 			grid.fixedNodes[:,grid.rightNodes]  .= fix
 			if ghostcell == true
           	   grid.fixedNodes[:,grid.rightNodes.-1]   .= fix
            end
 		else  error("Not recognized grid faces!!!")
 		end
 	end
end

# data["dirichlet_solid"] = [(1,"fix",(0,1,0)), # => fix  nodes of 'fix' group of solid 1 on Y dir
#                            (2,"fix",(0,1,0))] # => fix  nodes of 'fix' group of solid 2 on Y dir

function fix_Dirichlet_solid(solids,data)
 	if haskey(data, "dirichlet_solid") == false return end 

    dirichlet_solid = data["dirichlet_solid"]

 	for (s,bnd,fix) in dirichlet_solid
 		solid = solids[s]
        if haskey(solid.mesh.node_sets,bnd ) == false 
 		   Mesh.create_node_set_from_element_set!(solid.mesh, bnd)
 	    end
        ids  = collect(solid.mesh.node_sets[bnd])

 		solid.fixedNodes[:,ids] .= fix 		
 	end
end

# time-dependent Dirichlet boundary conditions on deformable solids
# data["time_dirichlet_solid"] = [(1,"fix",f)] # => fix  nodes of 'fix' group of solid 1 

function fix_Dirichlet_solid(solids,data,dt)
 	if haskey(data, "time_dirichlet_solid") == false return end 

    dirichlet_solid = data["time_dirichlet_solid"]

 	for (s,bnd,f) in dirichlet_solid
 		solid = solids[s]
        if haskey(solid.mesh.node_sets,bnd ) == false 
 		   Mesh.create_node_set_from_element_set!(solid.mesh, bnd)
 	    end
        ids  = collect(solid.mesh.node_sets[bnd])
        xx   = solid.pos
        XX   = solid.pos0
        du   = solid.dU
        vv   = solid.velocity

        # loop over the nodes of node set 'bnd', and apply boundary conditions defined by function 'f'
 		for id in ids
 			x0      = xx[id][1]
 			y0      = xx[id][2]
 			z0      = xx[id][3]
 			#println(sqrt(x0^2+y0^2))
 			(x,y,z,vx,vy,ux,uy)   = f(x0,y0,z0,dt)
 			xx[id]  = setindex(xx[id],x,1)
 			xx[id]  = setindex(xx[id],y,2)
 			xx[id]  = setindex(xx[id],z,3)
 			du[id]  = setindex(du[id],ux,1)
 			du[id]  = setindex(du[id],uy,2)
 			du[id]  = setindex(du[id],0.,3)
 			vv[id]  = setindex(vv[id],vx,1)
 			vv[id]  = setindex(vv[id],vy,2)
 			vv[id]  = setindex(vv[id],0.,3)
 		end	
 	end
end

################################################################################
# build_particle_neighbors: particle 'p' => neighbour particles
################################################################################

function build_particle_neighbors(solids,grid,rad)
	cellCount = (grid.nodeCountX-1)*(grid.nodeCountY-1)*(grid.nodeCountZ-1)
	for solid  in solids		
		xx                  = solid.pos		
		elems               = solid.elems				     
		feFuncs             = solid.N
		neighbours          = solid.neighbours


		centroids           = fill(zeros(3),solid.parCount)
		particle_in_cells   = [Set{Int64}() for _=1:cellCount] #fill(Set{Int64}(),cellCount)
	
	    # compute centroids of all elements
	  	@inbounds for ip = 1:solid.parCount
			elemNodes  =  @view elems[ip,:]  			
			coords     =  @view xx[elemNodes]
	    
			N          = feFuncs[ip]
            for i = 1 : length(N)
			  centroids[ip] += N[i] * coords[i]
		    end
		    c = get_cell(grid,centroids[ip])		    
		    push!(particle_in_cells[c],ip)
	  	end
	  	# build neighbors
	  	for ip = 1:solid.parCount
	  	  cells = get_cells(rad,centroids[ip],grid)
	  	  for cell in cells
	  	  	particles = particle_in_cells[cell]
	  	  	for p in particles
	  	  	   push!(neighbours[ip],p)
	  	    end
	  	  end
	  	end  
	end
end


export fix_Dirichlet_grid, fix_Dirichlet_solid, build_particle_neighbors

end