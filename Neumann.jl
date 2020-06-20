module Neumann

# This module implements functions for calculating external forces at nodes of the 
# FE meshes due to pressure.
# data is a dictionary: 
# data["pressure"] = [(1,"tag1",f),(2,"tag2",g)]
# means that there are 2 pressures: 
# f(t) applied on solid '1' (over element set 'tag1')
# g(t) applied on solid '2' (over element set 'tag2')

using StaticArrays
using Mesh

function compute_fext(solids,funcs_surface, normals_surface, weights_surface, gpCoords_surface,data,t)
 	if haskey(data, "pressure") == false return end 

 	pressure_data = data["pressure"]

    @inbounds for is = 1:length(pressure_data)
        solid_pressure_data = pressure_data[is]
        tag                 = solid_pressure_data[2]
        f                   = solid_pressure_data[3]
		# only deformable solids here
	  	solid  = solids[solid_pressure_data[1]]
	  	XX     = solid.pos0	  	  	  
	  	elems  = solid.mesh.elements
	  	basis  = solid.basis_S
	  	ftrac  = solid.ftrac
	  	for ip=1:solid.nodeCount 
	  		ftrac[ip]  = @SVector [0., 0., 0.]	  		
	  	end

        if haskey(solid.mesh.element_sets, tag)
	  	    surf_elems_ids = collect(solid.mesh.element_sets[tag])
	  	else
	  		error("No element set found")
	  	end
        # loop over  elements of the surface tag 'force'
	  	@inbounds for ip in surf_elems_ids
			elemNodes =  elems[ip]  
			coords    =  @view XX[elemNodes]
			#println(coords)
			getNormals!(funcs_surface, normals_surface, weights_surface , coords, gpCoords_surface, basis )
			# loop over Gauss points
			for ip=1:4                     
			    ww = weights_surface[ip]   
			    #println(normals_surface[:,ip])       
			    #println(ww)       
				for i = 1:length(elemNodes)
				   in  = elemNodes[i]; # index of node 'i'
				   ss  = funcs_surface[i,ip]*ww*f(t)#400*exp(-10000*t)
		           ftrac[in][1] -= normals_surface[1,ip]*ss		           
		           ftrac[in][2] -= normals_surface[2,ip]*ss
		           ftrac[in][3] -= normals_surface[3,ip]*ss
		        end
	        end
	   end
	end
end

export compute_fext	

end