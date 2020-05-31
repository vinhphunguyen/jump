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

module Mesh

using Printf
using LinearAlgebra
using StaticArrays
using DelimitedFiles

export lagrange_basis!, lagrange_basis_derivatives!
export Line2, Tri3, Quad4, Tet4, Hexa8

struct Line2 end
struct Tri3  end
struct Quad4 end
struct Hexa4 end
struct Tet4  end

function lagrange_basis!(N, dNdxi, type::Line2,coord)
	xi   = coord[1];
	N    = ([1-xi,1+xi]/2)';
	dNdxi= [-1;1]/2;
end

function lagrange_basis!(N, dNdxi, type::Tri3,coord)
	xi     = coord[1];
	eta    = coord[2];
	N      = [1-xi-eta;xi;eta];
	dNdxi  = [-1 -1;1. 0.;0. 1.]
end

function lagrange_basis!(N, dNdxi, type::Quad4,coord)
	xi     = coord[1];
	eta    = coord[2];
	N     .= 0.25*[ (1-xi)*(1-eta);
					        (1+xi)*(1-eta);
					        (1+xi)*(1+eta);
					        (1-xi)*(1+eta)];

  dNdxi.= 0.25 *[-(1-eta)   -(1-xi);
                   1-eta    -(1+xi);
                   1+eta      1+xi;
                 -(1+eta)     1-xi];
end

function lagrange_basis!(N, type::Quad4,xieta,coords)
	xi     = xieta[1];
	eta    = xieta[2];

	N      .=    SVector{4}( 0.25*(1-xi)*(1-eta),
					         0.25*(1+xi)*(1-eta),
					         0.25*(1+xi)*(1+eta),
					         0.25*(1-xi)*(1+eta) )
				         

    dN1dxi  = -0.25 * (1-eta)
    dN1deta = -0.25 * (1-xi)
    dN2dxi  =  0.25 * (1-eta)
    dN2deta = -0.25 * (1+xi)
    dN3dxi  =  0.25 * (1+eta)
    dN3deta =  0.25 * (1+xi)
    dN4dxi  = -0.25 * (1+eta)
    dN4deta =  0.25 * (1-xi)

    J11     = dN1dxi  * coords[1][1] + dN2dxi  * coords[2][1] + dN3dxi  * coords[3][1] + dN4dxi  * coords[4][1]
    J12     = dN1dxi  * coords[1][2] + dN2dxi  * coords[2][2] + dN3dxi  * coords[3][2] + dN4dxi  * coords[4][2]
    J21     = dN1deta * coords[1][1] + dN2deta * coords[2][1] + dN3deta * coords[3][1] + dN4deta * coords[4][1]
    J22     = dN1deta * coords[1][2] + dN2deta * coords[2][2] + dN3deta * coords[3][2] + dN4deta * coords[4][2]

    return J11*J22 - J12*J21

end

function lagrange_basis!(N, type::Hexa4,xieta,coords)
	xi     = xieta[1]
	eta    = xieta[2]
	zeta   = xieta[3]

	N      .=    SVector{8}( 0.125*(1-xi)*(1-eta)*(1-zeta),
					         0.125*(1+xi)*(1-eta)*(1-zeta),
					         0.125*(1+xi)*(1+eta)*(1-zeta),
					         0.125*(1-xi)*(1+eta)*(1-zeta),
					         0.125*(1-xi)*(1-eta)*(1+zeta),
					         0.125*(1+xi)*(1-eta)*(1+zeta),
					         0.125*(1+xi)*(1+eta)*(1+zeta),
					         0.125*(1-xi)*(1+eta)*(1+zeta)
					          )
				         

    dN1dxi   = -0.125*(1-eta)*(1-zeta)
    dN1deta  = -0.125*(1-xi)*(1-zeta)
    dN1dzeta = -0.125*(1-xi)*(1-eta)

    dN2dxi   =  0.125*(1-eta)*(1-zeta)
    dN2deta  = -0.125*(1+xi)*(1-zeta)
    dN2dzeta = -0.125*(1+xi)*(1-eta)

    dN3dxi   =  0.125*(1+eta)*(1-zeta)
    dN3deta  =  0.125*(1+xi)*(1-zeta)
    dN3dzeta = -0.125*(1+xi)*(1+eta)

    dN4dxi   = -0.125*(1+eta)*(1-zeta)
    dN4deta  =  0.125*(1-xi)*(1-zeta)
    dN4dzeta = -0.125*(1-xi)*(1+eta)

    dN5dxi   = -0.125*(1-eta)*(1+zeta)
    dN5deta  = -0.125*(1-xi)*(1+zeta)
    dN5dzeta = 0.125*(1-xi)*(1-eta)

    dN6dxi   = 0.125*(1-eta)*(1+zeta)
    dN6deta  = -0.125*(1+xi)*(1+zeta)
    dN6dzeta = 0.125*(1+xi)*(1-eta)

    dN7dxi   = 0.125*(1+eta)*(1+zeta)
    dN7deta  = 0.125*(1+xi)*(1+zeta)
    dN7dzeta = 0.125*(1+xi)*(1+eta)

    dN8dxi   = -0.125*(1+eta)*(1+zeta)
    dN8deta  =  0.125*(1-xi)*(1+zeta)
    dN8dzeta =  0.125*(1-xi)*(1+eta)

    J11 = dN1dxi * coords[1][1] + dN2dxi * coords[2][1] + dN3dxi  * coords[3][1] + dN4dxi * coords[4][1] +
          dN5dxi * coords[5][1] + dN6dxi * coords[6][1] + dN7dxi  * coords[7][1] + dN8dxi * coords[8][1]

    J12 = dN1dxi * coords[1][2] + dN2dxi * coords[2][2] + dN3dxi * coords[3][2] + dN4dxi * coords[4][2]+
          dN5dxi * coords[5][2] + dN6dxi * coords[6][2] + dN7dxi * coords[7][2] + dN8dxi * coords[8][2]

    J13 = dN1dxi * coords[1][3] + dN2dxi * coords[2][3] + dN3dxi * coords[3][3] + dN4dxi * coords[4][3]+
          dN5dxi * coords[5][3] + dN6dxi * coords[6][3] + dN7dxi * coords[7][3] + dN8dxi * coords[8][3]              
    J21 = dN1deta  * coords[1][1] + dN2deta  * coords[2][1] + dN3deta  * coords[3][1] + dN4deta * coords[4][1] +  dN5deta  * coords[5][1] + dN6deta  * coords[6][1] + dN7deta  * coords[7][1] + dN8deta * coords[8][1]

    J22 = dN1deta  * coords[1][2] + dN2deta  * coords[2][2] + dN3deta  * coords[3][2] + dN4deta * coords[4][2]+
          dN5deta  * coords[5][2] + dN6deta  * coords[6][2] + dN7deta  * coords[7][2] + dN8deta * coords[8][2]

    J23 = dN1deta  * coords[1][3] + dN2deta  * coords[2][3] + dN3deta  * coords[3][3] + dN4deta * coords[4][3]+
          dN5deta  * coords[5][3] + dN6deta  * coords[6][3] + dN7deta  * coords[7][3] + dN8deta * coords[8][3]    

    J31 = dN1dzeta * coords[1][1] + dN2dzeta * coords[2][1] + dN3dzeta * coords[3][1] + dN4dzeta * coords[4][1] + dN5dzeta * coords[5][1] + dN6dzeta * coords[6][1] + dN7dzeta * coords[7][1] + dN8dzeta * coords[8][1]

    J32 = dN1dzeta * coords[1][2] + dN2dzeta * coords[2][2] + dN3dzeta * coords[3][2] + dN4dzeta * coords[4][2]+
          dN5dzeta * coords[5][2] + dN6dzeta * coords[6][2] + dN7dzeta * coords[7][2] + dN8dzeta * coords[8][2]

    J33 = dN1dzeta * coords[1][3] + dN2dzeta * coords[2][3] + dN3dzeta * coords[3][3] + dN4dzeta * coords[4][3]+
          dN5dzeta * coords[5][3] + dN6dzeta * coords[6][3] + dN7dzeta * coords[7][3] + dN8dzeta * coords[8][3]      

    return det(SMatrix{3,3}(J11,J12,J13;J21,J22,J23;J31,J32,J33))

end

function lagrange_basis_derivatives!(dNdx, type::Quad4,xieta,coords)
	xi     = xieta[1];
	eta    = xieta[2];

    dN1dxi  = -0.25 * (1-eta)
    dN1deta = -0.25 * (1-xi)
    dN2dxi  =  0.25 * (1-eta)
    dN2deta = -0.25 * (1+xi)
    dN3dxi  =  0.25 * (1+eta)
    dN3deta =  0.25 * (1+xi)
    dN4dxi  = -0.25 * (1+eta)
    dN4deta =  0.25 * (1-xi)

    J11     = dN1dxi  * coords[1][1] + dN2dxi  * coords[2][1] + dN3dxi  * coords[3][1] + dN4dxi  * coords[4][1]
    J12     = dN1dxi  * coords[1][2] + dN2dxi  * coords[2][2] + dN3dxi  * coords[3][2] + dN4dxi  * coords[4][2]
    J21     = dN1deta * coords[1][1] + dN2deta * coords[2][1] + dN3deta * coords[3][1] + dN4deta * coords[4][1]
    J22     = dN1deta * coords[1][2] + dN2deta * coords[2][2] + dN3deta * coords[3][2] + dN4deta * coords[4][2]

    dNdx    .= inv(SMatrix{2,2}(J11,J21,J12,J22))*SMatrix{2,4}(dN1dxi,dN1deta,
    	                                                       dN2dxi,dN2deta,
    	                                                       dN3dxi,dN3deta,
    	                                                       dN4dxi,dN4deta) 
    return J11*J22 - J12*J21

end


function lagrange_basis_derivatives!(dNdx, type::Hexa4,xieta,coords)
	xi     = xieta[1]
	eta    = xieta[2]
	zeta   = xieta[3]

    dN1dxi   = -0.125*(1-eta)*(1-zeta)
    dN1deta  = -0.125*(1-xi)*(1-zeta)
    dN1dzeta = -0.125*(1-xi)*(1-eta)

    dN2dxi   =  0.125*(1-eta)*(1-zeta)
    dN2deta  = -0.125*(1+xi)*(1-zeta)
    dN2dzeta = -0.125*(1+xi)*(1-eta)

    dN3dxi   =  0.125*(1+eta)*(1-zeta)
    dN3deta  =  0.125*(1+xi)*(1-zeta)
    dN3dzeta = -0.125*(1+xi)*(1+eta)

    dN4dxi   = -0.125*(1+eta)*(1-zeta)
    dN4deta  =  0.125*(1-xi)*(1-zeta)
    dN4dzeta = -0.125*(1-xi)*(1+eta)

    dN5dxi   = -0.125*(1-eta)*(1+zeta)
    dN5deta  = -0.125*(1-xi)*(1+zeta)
    dN5dzeta = 0.125*(1-xi)*(1-eta)

    dN6dxi   = 0.125*(1-eta)*(1+zeta)
    dN6deta  = -0.125*(1+xi)*(1+zeta)
    dN6dzeta = 0.125*(1+xi)*(1-eta)

    dN7dxi   = 0.125*(1+eta)*(1+zeta)
    dN7deta  = 0.125*(1+xi)*(1+zeta)
    dN7dzeta = 0.125*(1+xi)*(1+eta)

    dN8dxi   = -0.125*(1+eta)*(1+zeta)
    dN8deta  =  0.125*(1-xi)*(1+zeta)
    dN8dzeta =  0.125*(1-xi)*(1+eta)

    J11 = dN1dxi * coords[1][1] + dN2dxi * coords[2][1] + dN3dxi  * coords[3][1] + dN4dxi * coords[4][1] +
          dN5dxi * coords[5][1] + dN6dxi * coords[6][1] + dN7dxi  * coords[7][1] + dN8dxi * coords[8][1]

    J12 = dN1dxi * coords[1][2] + dN2dxi * coords[2][2] + dN3dxi * coords[3][2] + dN4dxi * coords[4][2]+
          dN5dxi * coords[5][2] + dN6dxi * coords[6][2] + dN7dxi * coords[7][2] + dN8dxi * coords[8][2]

    J13 = dN1dxi * coords[1][3] + dN2dxi * coords[2][3] + dN3dxi * coords[3][3] + dN4dxi * coords[4][3]+
          dN5dxi * coords[5][3] + dN6dxi * coords[6][3] + dN7dxi * coords[7][3] + dN8dxi * coords[8][3]              
    J21 = dN1deta  * coords[1][1] + dN2deta  * coords[2][1] + dN3deta  * coords[3][1] + dN4deta * coords[4][1] +  dN5deta  * coords[5][1] + dN6deta  * coords[6][1] + dN7deta  * coords[7][1] + dN8deta * coords[8][1]

    J22 = dN1deta  * coords[1][2] + dN2deta  * coords[2][2] + dN3deta  * coords[3][2] + dN4deta * coords[4][2]+
          dN5deta  * coords[5][2] + dN6deta  * coords[6][2] + dN7deta  * coords[7][2] + dN8deta * coords[8][2]

    J23 = dN1deta  * coords[1][3] + dN2deta  * coords[2][3] + dN3deta  * coords[3][3] + dN4deta * coords[4][3]+
          dN5deta  * coords[5][3] + dN6deta  * coords[6][3] + dN7deta  * coords[7][3] + dN8deta * coords[8][3]    

    J31 = dN1dzeta * coords[1][1] + dN2dzeta * coords[2][1] + dN3dzeta * coords[3][1] + dN4dzeta * coords[4][1] + dN5dzeta * coords[5][1] + dN6dzeta * coords[6][1] + dN7dzeta * coords[7][1] + dN8dzeta * coords[8][1]

    J32 = dN1dzeta * coords[1][2] + dN2dzeta * coords[2][2] + dN3dzeta * coords[3][2] + dN4dzeta * coords[4][2]+
          dN5dzeta * coords[5][2] + dN6dzeta * coords[6][2] + dN7dzeta * coords[7][2] + dN8dzeta * coords[8][2]

    J33 = dN1dzeta * coords[1][3] + dN2dzeta * coords[2][3] + dN3dzeta * coords[3][3] + dN4dzeta * coords[4][3]+
          dN5dzeta * coords[5][3] + dN6dzeta * coords[6][3] + dN7dzeta * coords[7][3] + dN8dzeta * coords[8][3]      
    J  = SMatrix{3,3}(J11,J12,J13;J21,J22,J23;J31,J32,J33)      
    
    dNdx    .= inv(J)*SMatrix{3,8}(dN1dxi,dN1deta, dN1dzeta,
    	                           dN2dxi,dN2deta, dN2dzeta,
    	                           dN3dxi,dN3deta, dN3dzeta,
    	                           dN4dxi,dN4deta, dN4dzeta,
    	                           dN5dxi,dN5deta, dN5dzeta,
    	                           dN6dxi,dN6deta, dN6dzeta,
    	                           dN7dxi,dN7deta, dN7dzeta,
    	                           dN8dxi,dN8deta, dN8dzeta,
    	                           ) 

    return det(J)

end

function square_node_array(pt1,pt2,pt3,pt4,numnod_u,numnod_v)
    xi_pts  = LinRange(-1.,1.,numnod_u);
    eta_pts = LinRange(-1.,1.,numnod_v);
    x_pts   = [pt1[1],pt2[1],pt3[1],pt4[1]];
    y_pts   = [pt1[2],pt2[2],pt3[2],pt4[2]];

    N     = zeros(4)
	dNdxi = zeros(4,2)
	X     = zeros(2,numnod_u*numnod_v)
	for r=1:numnod_v
	  eta = eta_pts[r];
	  for c=1:numnod_u
	    xi = xi_pts[c];
	    # get interpolation basis at xi, eta
	    lagrange_basis!(N,dNdxi,Quad4(),[xi,eta])
	    X[:,(r-1)*numnod_u+c]=[dot(x_pts,N) dot(y_pts,N)]
		#println(N)
	  end
	end
	return X
end


function make_elem(num_u,num_v)
	nnx  = num_u+1
    nny  = num_v+1
	inc_u=1
	inc_v=nnx

	node_pattern = [1 2 nnx+2 nnx+1]
	inc          = zeros(Int64, 1,           size(node_pattern,2))
	element      = zeros(Int64, num_u*num_v, size(node_pattern,2))
    e = 1
	for row=1:num_v
	   for col=1:num_u
	      element[e,:]  = node_pattern+inc
	      inc         .+= inc_u;
	      e += 1;
	   end
	   inc .= row*inc_v
	end

	return element
end

function createMeshForRectangle(pt1,pt2,pt3,pt4,elemX,elemY)
        X       = square_node_array(pt1,pt2,pt3,pt4,elemX+1,elemY+1)
		element = make_elem(elemX,elemY)
		return (X,element)
end

# read Gmsh mesh file to get node coords and element connectivity
# In the mesh, there are multiple element types. For example, 2D mesh can have
# line elements (boundary) and solid elements.
# return (node coords,elem connectivity)
function loadGMSH(sFile::String)
		@printf("\nLoading .msh file: %s \n", sFile)
		hFile     = open(sFile)
		arrayLine = readlines(hFile)
		lengthCnt = length(arrayLine)

		# read the total number of nodes and elements from file
		nNodes::Int      = 0
		nElements::Int   = 0
		for indexLine = 1:lengthCnt
			if(occursin(arrayLine[indexLine], "\$Nodes"))
				nNodes = parse(Int, arrayLine[indexLine+1])
			end
			if(occursin(arrayLine[indexLine], "\$Elements"))
				nElements = parse(Int, arrayLine[indexLine+1])
			end
		end

		# reading node coords
		arrayNode_ID = Array{Int,1}(undef,nNodes)
		arrayNode_Coordinate = Array{Float64,2}(undef,3,nNodes)
		for indexLine = 1:lengthCnt
			if(occursin(arrayLine[indexLine], "\$Nodes"))
				@printf("	Reading %d nodes\n", nNodes)
				for indexNode = 1:nNodes
					arrayTemp = Array{Float64,1}(undef,4)
					arrayTemp = readdlm(IOBuffer(arrayLine[indexLine + 1 + indexNode])) # convert multi number string into array of numbers
					arrayNode_ID[indexNode] = arrayTemp[1]
					arrayNode_Coordinate[:,indexNode] = [arrayTemp[2]; arrayTemp[3]; arrayTemp[4]]
					# @printf("	Node_ID %d: Coordinates (%f,%f,%f)\n", arrayNode_ID[indexNode], arrayNode_Coordinate[indexNode,1], arrayNode_Coordinate[indexNode,2], arrayNode_Coordinate[indexNode,3])
				end
			end
		end

		# reading element data
		accepts = (3,5) # only 4 node quad and 8 node hexa elements
		arrayElement_ID = Array{Int,1}()
		inodes          = Array{Int,2}(undef,0,4) # for 4-node tetrahedron/quad only
		for indexLine = 1:lengthCnt
			if(occursin(arrayLine[indexLine], "\$Elements"))
				@printf("	Reading %d elements\n", nElements)
				for indexElement = 1:nElements
					# from gmsh document -> elm-number elm-type number-of-tags < tag > ... node-number-list
					arrayTemp = readdlm(IOBuffer(arrayLine[indexLine + 1 + indexElement]),Int16) # convert multi number string into array of numbers
					elemType  = arrayTemp[2]
					if !(elemType in accepts) # element type = 4-node quad
						continue
					end
					if(elemType == 3) # element type = 4-node quad
						#println(arrayTemp)
						push!(arrayElement_ID,arrayTemp[1])
						nTags           = arrayTemp[3]
						newCorners      = [arrayTemp[3+nTags+1] arrayTemp[3+nTags+2] arrayTemp[3+nTags+3] arrayTemp[3+nTags+4]]
                        #println(newCorners)
						inodes = [inodes;newCorners]
						#println(inodes)
					end
				end
			end
		end
		#println(arrayElement_ID)
        #println(inodes)

		@printf("	nNodes: %d, nElements: %d\n", nNodes, nElements)

		return (arrayNode_Coordinate,inodes)
	end

export createMeshForRectangle, loadGMSH

end
