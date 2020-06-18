module Fem

using StaticArrays
using  WriteVTK
using Printf

using Material
using Mesh
using Grid
using Util
using AbaqusReader

###########################################################
# Struct for 2D solid: a set of data for 2D particles
# Createa solid by: solid = Solid2D(coords,mat)
# Createa solid by: solid = Solid2D(coords,mat)
###########################################################
struct FEM2D
	mass                :: Vector{Float64}
	volume              :: Vector{Float64}
	centerX             :: Vector{Float64}             # centroids of all elements in original configuration

	pos0                :: Vector{SVector{2,Float64}}  # node coords initial
	pos                 :: Vector{SVector{2,Float64}}  # position
	velocity            :: Vector{SVector{2,Float64}}  # velocity
	dU                  :: Vector{SVector{2,Float64}}  # incremental displacements
	fint                :: Vector{SVector{2,Float64}}  # internal forces at FE nodes
	fbody               :: Vector{SVector{2,Float64}}  # external forces  due to gravity at FE nodes
	ftrac               :: Vector{SVector{2,Float64}}  # external forces  due to traction/pressure at FE nodes

	deformationGradient :: Vector{SMatrix{2,2,Float64,4}}  # F, 2x2 matrix
	strain              :: Vector{SMatrix{2,2,Float64,4}}  # stress, 2x2 matrix
	stress              :: Vector{MMatrix{2,2,Float64,4}}  # strain
	

	parCount            :: Int64              # same as element count in the mesh
	nodeCount           :: Int64              # number of nodes in the mesh

	elems               :: Array{Int64,2}
	mesh                :: FEMesh
	

    fixedNodesX         :: Vector{Int64}
	fixedNodesY         :: Vector{Int64}

	basis               :: FiniteElement     # finite element basis function
	vtk_cell            :: VTKCellType
	#boundary_nodes      :: Set{Int64}

	# particles from a mesh
	function FEM2D(fileName)
		mesh      = read_GMSH(fileName)
	    nodeCount = length(mesh.nodes)                       # node count
		parCount  = length(mesh.element_sets["All"])         # element count
		Identity = SMatrix{2,2}(1, 0, 0, 1)
		F        = fill(Identity,parCount)
		strain   = fill(zeros(2,2),parCount)
		stress   = fill(zeros(2,2),parCount)
		vol      = fill(0.,parCount)
		m        = fill(0.,nodeCount)
		x        = fill(zeros(2),nodeCount)
		nodesX   = fill(zeros(2),nodeCount)
		centerX  = fill(0.,parCount)
		elems    = Array{Array{Int64,1},1}(undef,0)   # element nodes
		velo     = fill(zeros(2),nodeCount)
	
		#println(size(nodes,2))
		println(parCount)

        # convert from mesh.nodes to our traditional data structure
		for i=1:nodeCount			
			nodesX[i] = @view mesh.nodes[i][1:2]
		end

        # convert from mesh.elements to our traditional data structure
        volumetricElems = collect(mesh.element_sets["All"])
		for i=1:parCount			
			push!(elems,mesh.elements[volumetricElems[i]])
		end

   
		nnodePerElem = length(elems[1])


		if      nnodePerElem == 3 
			basis = Tri3()  
			vtk_cell = VTKCellTypes.VTK_TRIANGLE
		elseif  nnodePerElem == 4 
			basis = Quad4() 
			vtk_cell = VTKCellTypes.VTK_QUAD
		else
			@printf("Number of nodes per element: %d \n", nnodePerElem)
			error("Unsupported element types: only Tri3 and Quad4 now!\n")
		end

		# if nodePerElem == 4 
		# 	meshBasis = Tet4()  
		# 	wgt       = .166666667
		# 	weights   = [0.166666667]   # this is so dangerous!!! all books ay weight =1
		# 	gpCoords  = [0.25,0.25,0.25]
	 #    end

	    fixX       = fill(0,nodeCount)
		fixY       = fill(0,nodeCount)

		if (haskey(mesh.element_sets,"boundary")) Mesh.create_node_set_from_element_set!(mesh, "boundary") end
		

		new(m,vol,centerX,nodesX,copy(nodesX),velo,copy(velo),copy(velo),copy(velo),copy(velo),F,
			   strain,stress,parCount,nodeCount,vcat(map(x->x', elems)...),mesh, fixX, fixY, basis, vtk_cell)
	end
end

###########################################################
# Struct for 3D solid: a set of data for 2D particles
# Createa solid by: solid = FEM3D("file.msh")  <= GMsh
# Createa solid by: solid = FEM3D("file.inp")  <= Abaqus
###########################################################
struct FEM3D
	mass                :: Vector{Float64}
	volume              :: Vector{Float64}

	pos0                :: Vector{SVector{3,Float64}}  # node coords initial
	pos                 :: Vector{SVector{3,Float64}}  # position
	velocity            :: Vector{SVector{3,Float64}}  # velocity
	dU                  :: Vector{SVector{3,Float64}}  # incremental displacements
	fint                :: Vector{SVector{3,Float64}}  # internal forces at FE nodes
	fbody               :: Vector{SVector{3,Float64}}  # external forces  due to gravity at FE nodes
	ftrac               :: Vector{MVector{3,Float64}}  # external forces  due to traction/pressure at FE nodes

	deformationGradient :: Vector{SMatrix{3,3,Float64,9}}  # F, 2x2 matrix
	strain              :: Vector{SMatrix{3,3,Float64,9}}  # stress, 2x2 matrix
	stress              :: Vector{MMatrix{3,3,Float64,9}}  # strain
	

	parCount            :: Int64              # same as element count in the mesh
	nodeCount           :: Int64              # number of nodes in the mesh	

	elems               :: Array{Int64,2}
	mesh                :: FEMesh
	fixedNodesX         :: Vector{Int64}
	fixedNodesY         :: Vector{Int64}
	fixedNodesZ         :: Vector{Int64}

	basis               :: FiniteElement
	vtk_cell            :: VTKCellType

    detJ                :: Vector{Float64} 
    dNdx                :: Vector{Array{Float64,2}}  # derivatives of all shape functions
    N                   :: Vector{Vector{Float64}}  # shape functions of all elements

	# particles from a mesh
	function FEM3D(fileName)

		ext       = fileName[findlast(isequal('.'),fileName) : end]

		if     ext == ".msh" 
			mesh      = read_GMSH(fileName) 
			nodeCount = length(mesh.nodes)                       # node count
		    parCount  = length(mesh.element_sets["All"])         # element count
		    nodesX    = fill(zeros(3),nodeCount)
		    elems     = Array{Array{Int64,1},1}(undef,0)   # element nodes
	        # convert from mesh.nodes to our traditional data structure
			for i=1:nodeCount			
				nodesX[i] = mesh.nodes[i]
			end

	        # convert from mesh.elements to our traditional data structure
	        volumetricElems = collect(mesh.element_sets["All"])
			for i=1:parCount			
				push!(elems,mesh.elements[volumetricElems[i]])
			end
		elseif ext == ".inp"
			mesh_abaqus     = abaqus_read_mesh(fileName)
			nodeCount       = length(mesh_abaqus["nodes"])
			elem_sets       = mesh_abaqus["element_sets"]
			volumetricElems = elem_sets["ALL"]
			parCount        = length(volumetricElems)
			elements        = mesh_abaqus["elements"]
            nodesX          = fill(zeros(3),nodeCount)
			elems           = Array{Array{Int64,1},1}(undef,0)   # element nodes
			nodes           = mesh_abaqus["nodes"]
			for i=1:nodeCount			
				nodesX[i] = nodes[i]
			end

			for i=1:parCount			
				push!(elems,elements[volumetricElems[i]])
			end

			mesh = FEMesh()
		end


		Identity  = SMatrix{3,3}(1, 0, 0,0, 1, 0,0, 0, 1)
		F         = fill(Identity,parCount)
		strain    = fill(zeros(3,3),parCount)
		stress    = fill(zeros(3,3),parCount)
		vol       = fill(0.,parCount)
		detJ      = fill(0.,parCount)
		m         = fill(0.,nodeCount)

		
		x         = fill(zeros(3),nodeCount)
		
		
		
		velo      = fill(zeros(3),nodeCount)
		#println(size(nodes,2))
		#println(parCount)


		fixX       = fill(0,nodeCount)
		fixY       = fill(0,nodeCount)
		fixZ       = fill(0,nodeCount)


		nnodePerElem = length(elems[1])

		if nnodePerElem == 4 
			basis     = Tet4()  
			vtk_cell  = VTKCellTypes.VTK_TETRA
			dNdx      = fill(zeros(3,4),parCount)
            N         = fill(zeros(4),parCount)
		end
		if nnodePerElem == 8 
			basis = Hexa8() 
			vtk_cell  = VTKCellTypes.VTK_HEXAHEDRON
			dNdx      = fill(zeros(3,8),parCount)
            N         = fill(zeros(8),parCount)
		end



		new(m,vol,nodesX,copy(nodesX),velo,copy(velo),copy(velo),copy(velo),copy(velo),F,
			   strain,stress,parCount,nodeCount,vcat(map(x->x', elems)...), mesh, fixX,fixY,fixZ,basis,vtk_cell,detJ,dNdx,N)
	end
end

###########################################################
# Struct for axi-symmetric solid: a set of data for 2D particles
# Createa solid by: solid = FEMAxis("file.msh")  <= GMsh
# Createa solid by: solid = FEMAxis("file.inp")  <= Abaqus
###########################################################
struct FEMAxis
	mass                :: Vector{Float64}
	volume              :: Vector{Float64}

	pos0                :: Vector{SVector{2,Float64}}  # node coords initial
	pos                 :: Vector{SVector{2,Float64}}  # position
	velocity            :: Vector{SVector{2,Float64}}  # velocity
	dU                  :: Vector{SVector{2,Float64}}  # incremental displacements
	fint                :: Vector{SVector{2,Float64}}  # internal forces at FE nodes
	fbody               :: Vector{SVector{2,Float64}}  # external forces  due to gravity at FE nodes
	ftrac               :: Vector{MVector{2,Float64}}  # external forces  due to traction/pressure at FE nodes

	deformationGradient :: Vector{SMatrix{3,3,Float64,9}}  # F, 2x2 matrix
	strain              :: Vector{SMatrix{3,3,Float64,9}}  # stress, 2x2 matrix
	stress              :: Vector{MMatrix{3,3,Float64,9}}  # strain
	

	parCount            :: Int64              # same as element count in the mesh
	nodeCount           :: Int64              # number of nodes in the mesh	

	elems               :: Array{Int64,2}
	mesh                :: FEMesh
	fixedNodesX         :: Vector{Int64}
	fixedNodesY         :: Vector{Int64}

	basis               :: FiniteElement
	vtk_cell            :: VTKCellType

    detJ                :: Vector{Float64} 
    dNdx                :: Vector{Array{Float64,2}}  # derivatives of all shape functions
    N                   :: Vector{Vector{Float64}}  # shape functions of all elements

	# particles from a mesh
	function FEMAxis(fileName)

		ext       = fileName[findlast(isequal('.'),fileName) : end]

		if     ext == ".msh" 
			mesh      = read_GMSH(fileName) 
			nodeCount = length(mesh.nodes)                       # node count
		    parCount  = length(mesh.element_sets["All"])         # element count
		    nodesX    = fill(zeros(2),nodeCount)
		    elems     = Array{Array{Int64,1},1}(undef,0)   # element nodes
	        # convert from mesh.nodes to our traditional data structure
			for i=1:nodeCount			
				nodesX[i] = @view mesh.nodes[i][1:2]
			end

	        # convert from mesh.elements to our traditional data structure
	        volumetricElems = collect(mesh.element_sets["All"])
			for i=1:parCount			
				push!(elems,mesh.elements[volumetricElems[i]])
			end
		elseif ext == ".inp"
			mesh_abaqus     = abaqus_read_mesh(fileName)
			nodeCount       = length(mesh_abaqus["nodes"])
			elem_sets       = mesh_abaqus["element_sets"]
			volumetricElems = elem_sets["ALL"]
			parCount        = length(volumetricElems)
			elements        = mesh_abaqus["elements"]
            nodesX          = fill(zeros(2),nodeCount)
			elems           = Array{Array{Int64,1},1}(undef,0)   # element nodes
			nodes           = mesh_abaqus["nodes"]
			for i=1:nodeCount			
				nodesX[i] = nodes[i]
			end

			for i=1:parCount			
				push!(elems,elements[volumetricElems[i]])
			end

			mesh = FEMesh()
		end


		Identity  = SMatrix{3,3}(1, 0, 0,0, 1, 0,0, 0, 1)
		F         = fill(Identity,parCount)
		strain    = fill(zeros(3,3),parCount)
		stress    = fill(zeros(3,3),parCount)
		vol       = fill(0.,parCount)
		detJ      = fill(0.,parCount)
		m         = fill(0.,nodeCount)

		
		x         = fill(zeros(2),nodeCount)						
		velo      = fill(zeros(2),nodeCount)
		#println(size(nodes,2))
		#println(parCount)


		fixX       = fill(0,nodeCount)
		fixY       = fill(0,nodeCount)


		nnodePerElem = length(elems[1])

		if nnodePerElem == 3 
			basis     = Tri3()  
			vtk_cell  = VTKCellTypes.VTK_TRIANGLE
			dNdx      = fill(zeros(2,3),parCount)
            N         = fill(zeros(3),parCount)
		end
		if nnodePerElem == 4
			basis     = Quad4() 
			vtk_cell  = VTKCellTypes.VTK_QUAD
			dNdx      = fill(zeros(2,4),parCount)
            N         = fill(zeros(4),parCount)
		end

        #println(elems)
		new(m,vol,nodesX,copy(nodesX),velo,copy(velo),copy(velo),copy(velo),copy(velo),F,
			   strain,stress,parCount,nodeCount,vcat(map(x->x', elems)...), mesh, fixX,fixY,basis,vtk_cell,detJ,dNdx,N)
	end
end

function move(solid,dx)
   x  = solid.pos
   x0 = solid.pos0
   @inbounds for p = 1 : solid.nodeCount
	  x[p]  += dx
	  x0[p] += dx
   end
end

function rotate(solid,alpha)
   xx = solid.pos
   x0 = solid.pos0
   beta = deg2rad(alpha)
   @inbounds for p = 1 : solid.nodeCount
	  xp = xx[p][1]
	  yp = xx[p][2]
	  xq = xp*cos(beta)-yp*sin(beta)
	  yq = yp*cos(beta)+xp*sin(beta)
	  xx[p] = setindex(xx[p],xq,1)
	  xx[p] = setindex(xx[p],yq,2)
	  x0[p] = setindex(x0[p],xq,1)
	  x0[p] = setindex(x0[p],yq,2)
   end
end

function assign_velocity(solid,v0)
       velo = solid.velocity
       @inbounds for p = 1 : solid.nodeCount
          velo[p] = v0
       end
end

function fixYNodes(solid, tag)
  Mesh.create_node_set_from_element_set!(solid.mesh, tag)
  ids  = collect(solid.mesh.node_sets[tag])
  solid.fixedNodesY[ids] .= 1
end

function fixNodes(solid, tag)
  Mesh.create_node_set_from_element_set!(solid.mesh, tag)
  ids  = collect(solid.mesh.node_sets[tag])
  solid.fixedNodesX[ids] .= 1
  solid.fixedNodesY[ids] .= 1
end

function initializeBasis(solid::FEMAxis,density)
	nodePerElem = size(solid.elems,2)

	if (nodePerElem == 3)  		
		weights   = 0.5   # this is so dangerous!!! all books ay weight =1
		gpCoords  = [0.3333333333333,0.3333333333333]
		funcs     = zeros(3)#@SVector [0,0,0,0]
		ders      = zeros(2,3)
    end

    if (nodePerElem == 4)  	    
		weights   = 4.0   # this is so dangerous!!! all books ay weight =1
		gpCoords  = [0.0,0.0]
		funcs     = zeros(4)#@SVector [0,0,0,0]
		ders      = zeros(2,4)
	end

	xx     = solid.pos
	mm     = solid.mass  # to be updated here
	elems  = solid.elems		
	vol    = solid.volume
	meshBasis = solid.basis
	N      = solid.N
	detJ   = solid.detJ
	dNdx   = solid.dNdx

  	@inbounds for ip = 1:solid.parCount
		elemNodes  =  @view elems[ip,:]  
		elemNodes0 =        elems[ip,:]  
		coords     =  @view xx[elemNodes]
    
		J          = lagrange_basis_derivatives!(funcs, ders, meshBasis, gpCoords, coords)
		if J < 0.
		    if     nodePerElem == 4
		  	    elemNodes[2] = elemNodes0[4]			  	
		  	    elemNodes[4] = elemNodes0[2]
		    elseif nodePerElem == 3
                elemNodes[2] = elemNodes0[3]			  	
		  	    elemNodes[3] = elemNodes0[2]
		    end
		    coords     =  @view xx[elemNodes]
		    J          = lagrange_basis_derivatives!(funcs, ders, meshBasis, gpCoords, coords)
		end

        # only Quad4, but we do not want to use Tri3, do we?
		r         = funcs[1] * coords[1][1] + funcs[2] * coords[2][1] + funcs[3] * coords[3][1] + funcs[4] * coords[4][1]

        www       = r*J*weights
		detJ[ip]  = www
		N[ip]     = copy(funcs/r)  # be careful with this, as funcs,ders change from element to element
		dNdx[ip]  = copy(ders)
		
		vol[ip]   = www
		for i=1:length(elemNodes)
		  	id      = elemNodes[i]
		  	mm[id]  += density*funcs[i]*www			  	
		end
  	end
end	

function initializeBasis(solid::FEM3D,density)
	nodePerElem = size(solid.elems,2)

	if (nodePerElem == 4)  		
		weights   =  0.166666667  # this is so dangerous!!! all books ay weight =1
		gpCoords  = [0.25,0.25,0.25]
		funcs     = zeros(4)
        ders      = zeros(3,4)
    end

    if (nodePerElem == 8)  	    
        weights   = 8.
        gpCoords  = [0.,0.,0.]
        funcs     = zeros(8)
        ders      = zeros(3,8)
	end

	xx     = solid.pos
	mm     = solid.mass  # to be updated here
	elems  = solid.elems		
	vol    = solid.volume
	meshBasis = solid.basis
	N      = solid.N
	detJ   = solid.detJ
	dNdx   = solid.dNdx

  	@inbounds for ip = 1:solid.parCount
		elemNodes  =  @view elems[ip,:]  
		elemNodes0 =        elems[ip,:]  
		coords     =  @view xx[elemNodes]
    
		J          = lagrange_basis_derivatives!(funcs, ders, meshBasis, gpCoords, coords)
		if J < 0.
		  	#elemNodes[2] = elemNodes0[4]			  	
		  	#elemNodes[4] = elemNodes0[2]
		  	error("Mesh with negative Jacobian!")
		end

        www       = J*weights
		detJ[ip]  = www
		N[ip]     = copy(funcs)  # be careful with this, as funcs,ders change from element to element
		dNdx[ip]  = copy(ders)
		
		vol[ip]   = www
		for i=1:length(elemNodes)
		  	id      = elemNodes[i]
		  	mm[id]  += density*funcs[i]*www			  	
		end
  	end
end	

export FEM2D, FEM3D, FEMAxis, move, assign_velocity, fixYNodes, fixNodes, rotate, initializeBasis


end