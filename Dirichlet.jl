module Dirichlet

# This module implements functions for applying Dirichlet boundary conditions on the Eulerian grid
# data["dirichlet_grid"] = [("bottom",(0,1,0)),      # => fix bottom nodes on Y dir
#                           ("top",(1,0,0)),      # => fix top nodes on X dir
#                          ]

using StaticArrays
using Mesh

function fix_Dirichlet_grid(grid,data)
 	if haskey(data, "dirichlet_grid") == false return end 

    dirichlet_grid = data["dirichlet_grid"]

 	for (bnd,fix) in dirichlet_grid

 		if     bnd == "bottom" grid.fixedNodes[:,grid.bottomNodes] .= fix
 		elseif bnd == "top"    grid.fixedNodes[:,grid.topNodes]    .= fix
 		elseif bnd == "left"   grid.fixedNodes[:,grid.leftNodes]   .= fix
 		elseif bnd == "right"  grid.fixedNodes[:,grid.rightNodes]  .= fix
 		elseif bnd == "front"  grid.fixedNodes[:,grid.frontNodes]  .= fix
 		elseif bnd == "back"   grid.fixedNodes[:,grid.bottomNodes] .= fix
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

export fix_Dirichlet_grid, fix_Dirichlet_solid

end