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

module Fem

using  StaticArrays
using  Statistics
using  WriteVTK
using  Printf
using  WriteVTK

using Material
using Mesh
using Grid
using Util
using AbaqusReader

###########################################################
# Struct for 2D solid: a set of data for 2D FE particles
# Createa solid by: solid = FEM2D("file.msh")
# Createa solid by: solid = FEM2D("file.inp")
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
    stress              :: Vector{SMatrix{2,2,Float64,4}}  # strain
    

    parCount            :: Int64              # same as element count in the mesh
    nodeCount           :: Int64              # number of nodes in the mesh

    elems               :: Array{Int64,2}
    mesh                :: FEMesh
    
    fixedNodes          :: Array{Int64,2}

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

        fixedNodes   = Array{Int64}(undef, 2, nodeCount)
        fixedNodes  .= 0

	if (haskey(mesh.element_sets,"boundary")) Mesh.create_node_set_from_element_set!(mesh, "boundary") end
	

	new(m,vol,centerX,nodesX,copy(nodesX),velo,copy(velo),copy(velo),copy(velo),copy(velo),F,
	    strain,stress,parCount,nodeCount,vcat(map(x->x', elems)...),mesh, fixedNodes, basis, vtk_cell)
    end
end

###########################################################
# Struct for 2D rigid solid: a set of data for 2D FE particles
# Createa solid by: solid = Rigid2D("file.msh")
# Createa solid by: solid = Rigid2D("file.inp")
###########################################################
struct Rigid2D
    pos0                :: Vector{SVector{2,Float64}}  # position
    pos                 :: Vector{SVector{2,Float64}}  # position    
    centroids           :: Vector{SVector{2,Float64}}  # centroids of boundary elems
    normals             :: Vector{SVector{2,Float64}}  # normal vectors
    velocity            :: Vector{SVector{2,Float64}}  # velocity
    stress              :: Vector{MMatrix{2,2,Float64,4}}  # for compability with FEM2D/3D (for output)

    parCount            :: Int64              # same as element count in the mesh
    nodeCount           :: Int64              # number of nodes in the mesh

    elems               :: Array{Int64,2}
    surface_elems       :: Array{Int64,2}
    mesh                :: FEMesh
    basis               :: FiniteElement     # finite element basis function
    vtk_cell            :: VTKCellType
    #boundary_nodes      :: Set{Int64}

    # particles from a mesh
    function Rigid2D(fileName)
	mesh          = read_GMSH(fileName)
	boundaryElems = collect(mesh.element_sets["boundary"])
	solidElems    = collect(mesh.element_sets["All"])
	parCount      = length(solidElems)         # element count
	elems         = Array{Array{Int64,1},1}(undef,0)   # element nodes
	surface_elems = Array{Array{Int64,1},1}(undef,0)   # element nodes
	for i=1:parCount			
	    push!(elems,mesh.elements[solidElems[i]])
	end
	for i=1:length(boundaryElems)			
	    push!(surface_elems,mesh.elements[boundaryElems[i]])
	end

	nodeCount       = length(mesh.nodes)                       # node count
	nodes           = fill(zeros(2),nodeCount)
	centroids       = fill(zeros(2),length(boundaryElems)	)
	normals         = fill(zeros(2),length(boundaryElems)	)
	
	# convert from mesh.nodes to our traditional data structure
	for i=1:nodeCount			
	    nodes[i] = @view mesh.nodes[i][1:2]
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

	if (haskey(mesh.element_sets,"boundary")) Mesh.create_node_set_from_element_set!(mesh, "boundary") end

	stress   = fill(zeros(2,2),parCount)
	velo     = fill(zeros(2),nodeCount)
	
	new(nodes, copy(nodes), centroids, normals, velo, stress, parCount,nodeCount,
	    vcat(map(x->x', elems)...),vcat(map(x->x', surface_elems)...), mesh, basis, vtk_cell)
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
    ftrac               :: Vector{SVector{3,Float64}}  # external forces  due to traction/pressure at FE nodes

    deformationGradient :: Vector{SMatrix{3,3,Float64,9}}  # F, 2x2 matrix
    strain              :: Vector{SMatrix{3,3,Float64,9}}  # strain, 3x3 matrix
    stress              :: Vector{SMatrix{3,3,Float64,9}}  # stress


    parCount            :: Int64              # same as element count in the mesh
    nodeCount           :: Int64              # number of nodes in the mesh	
    surfCount           :: Int64              # number of centroid

    elems               :: Array{Int64,2}
    mesh                :: FEMesh
    fixedNodes          :: Array{Int64,2}

    basis               :: FiniteElement
    basis_S             :: FiniteElement
    vtk_cell            :: VTKCellType

    detJ                :: Vector{Float64} 
    dNdx                :: Vector{Array{Float64,2}}  # derivatives of all shape functions
    N                   :: Vector{Vector{Float64}}  # shape functions of all elements

    reaction_forces     :: MVector{3,Float64}
    neighbours          :: Vector{Set{Int64}}       # neighbours element for a give elem
    elem_mass           :: Vector{Float64}

    # for rigid bodies

    #centroids           :: Vector{SVector{3,Float64}}  # centroids of boundary elems
    normals             :: Vector{SVector{3,Float64}}  # normal vectors
    surface_elems       :: Array{Int64,2}

    # VTK output related

    cells               :: Vector{MeshCell}

    # particles from a mesh
    function FEM3D(fileName)

	ext       = fileName[findlast(isequal('.'),fileName) : end]

	if     ext == ".msh" 
	    mesh      = read_GMSH(fileName) 
	    nodeCount = length(mesh.nodes)                       # node count
	    parCount  = length(mesh.element_sets["All"])         # element count
	    nodesX    = fill(zeros(3),nodeCount)
	    elems     = Array{Array{Int64,1},1}(undef,0)   # element nodes
	    surface_elems= Array{Array{Int64,1},1}(undef,0)   # element nodes
	    # convert from mesh.nodes to our traditional data structure
	    for i=1:nodeCount			
		nodesX[i] = mesh.nodes[i]
	    end

            neighbours      = Vector{Set{Int64}}(undef,parCount)
	    # convert from mesh.elements to our traditional data structure
	    volumetricElems = collect(mesh.element_sets["All"])
	    
	    for i=1:parCount			
		push!(elems,mesh.elements[volumetricElems[i]])		
		push!(neighbours,Set{Int64}())	
	    end
	    
	    
	    if (haskey(mesh.element_sets,"boundary")) 
		boundaryElems   = collect(mesh.element_sets["boundary"])
		for i=1:length(boundaryElems)			
		    push!(surface_elems,mesh.elements[boundaryElems[i]])
		end

		Mesh.create_node_set_from_element_set!(mesh, "boundary") 
		#centroids       = fill(zeros(3),length(mesh.node_sets["boundary"])	)
		surfCount       = length(mesh.node_sets["boundary"])
            else
		centroids       = fill(zeros(3),0)
		surface_elems   = [[1,1],[1,1]]
		surfCount       = 1
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
	    
	    neighbours      = Vector{Set{Int64}}(undef,parCount)

	    for i=1:nodeCount			
		nodesX[i] = nodes[i]
	    end

	    for i=1:parCount			
		push!(elems,elements[volumetricElems[i]])
		push!(neighbours,Set{Int64}())
	    end

	    mesh = FEMesh()
	    mesh.elements = elements
	    # convert from mesh_abaqus to mesh (use slightly different data structure)
	    #mesh.node_sets = mesh_abaqus["node_sets"]
	    for (key,val) in mesh_abaqus["node_sets"]
		node_ids = Set{Int}()
		for elid in val
		    push!(node_ids, elid)
		end
		mesh.node_sets[key] = node_ids
	    end
	    for (key,val) in mesh_abaqus["element_sets"]
		node_ids = Set{Int}()
		for elid in val
		    push!(node_ids, elid)
		end
		mesh.element_sets[key] = node_ids
	    end

            if (haskey(elem_sets,"boundary")) 
            	boundaryElems   = elem_sets["boundary"]
            	#centroids       = fill(zeros(3),length(mesh.node_sets["boundary"])	)
            	surfCount       = length(mesh.node_sets["boundary"])
            	surface_elems   = Array{Array{Int64,1},1}(undef,0)   # element nodes	
		for i=1:length(boundaryElems)			
		    push!(surface_elems,elements[boundaryElems[i]])
		end
	    else
		#centroids       = fill(zeros(3),0)
		surface_elems   = [[1,1],[1,1]]
		surfCount       = 1
	    end
	end


	Identity  = SMatrix{3,3}(1, 0, 0,0, 1, 0,0, 0, 1)
	F         = fill(Identity,parCount)
	strain    = fill(zeros(3,3),parCount)
	stress    = fill(zeros(3,3),parCount)
	vol       = fill(0.,parCount)
	detJ      = fill(0.,parCount)
	m         = fill(0.,nodeCount)
	em        = fill(0.,parCount)

	x         = fill(zeros(3),nodeCount)
	velo      = fill(zeros(3),nodeCount)
	#println(size(nodes,2))
	#println(parCount)

        fixedNodes   = Array{Int64}(undef, 3, nodeCount)
        fixedNodes  .= 0


	nnodePerElem = length(elems[1])

	if nnodePerElem == 4 
	    basis     = Tet4()  
	    basis_S   = Tri3()  
	    vtk_cell  = VTKCellTypes.VTK_TETRA
	    dNdx      = fill(zeros(3,4),parCount)
            N         = fill(zeros(4),parCount)
	end
	if nnodePerElem == 8 
	    basis     = Hexa8() 
	    basis_S   = Quad4() 
	    vtk_cell  = VTKCellTypes.VTK_HEXAHEDRON
	    dNdx      = fill(zeros(3,8),parCount)
            N         = fill(zeros(8),parCount)
	end

	cells = Vector{MeshCell}(undef,parCount)

	for e=1:parCount
	    inds = elems[e] 
	    c    = MeshCell(vtk_cell, inds)
            cells[e] = c
        end

	new(m,vol,nodesX,copy(nodesX),velo,copy(velo),copy(velo),copy(velo),copy(velo),F,
	    strain,stress,parCount,nodeCount,surfCount,vcat(map(x->x', elems)...), mesh, 
	    fixedNodes,basis,basis_S,vtk_cell,detJ,dNdx,N,zeros(3),neighbours,em,
	    copy(nodesX),vcat(map(x->x', surface_elems)...), cells )
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
    me     = solid.elem_mass

    @inbounds for ip = 1:solid.parCount
	elemNodes  =  @view elems[ip,:]  
	elemNodes0 =        elems[ip,:]  
	coords     =  @view xx[elemNodes]
        
	J          = lagrange_basis_derivatives!(funcs, ders, meshBasis, gpCoords, coords)
	if J < 0.
	    #elemNodes[2] = elemNodes0[4]			  	
	    #elemNodes[4] = elemNodes0[2]
	    println(meshBasis)
	    error("Mesh with negative Jacobian!")
	end

        www       = J*weights
	detJ[ip]  = www
	N[ip]     = copy(funcs)  # be careful with this, as funcs,ders change from element to element
	dNdx[ip]  = copy(ders)
	
	vol[ip]   = www
	me[ip]    = www * density
	for i=1:length(elemNodes)
	    id      = elemNodes[i]
	    mm[id]  += density*funcs[i]*www			  	
	end
    end
end	


function compute_normals!(solid::Rigid2D)
    elems   = solid.surface_elems
    xx      = solid.pos
    normals = solid.normals
    centroids = solid.centroids
    normal  = zeros(2)
    @inbounds for ip = 1:size(elems,1)
	elemNodes  =  @view elems[ip,:]  		
	coords     =  @view xx[elemNodes]
	getNormals!(normal,coords,Line2())
	normals[ip]   += normal
	centroids[ip] += 0.5*@SVector[coords[1][1]+coords[2][1],coords[1][2]+coords[2][2]]
    end		
end

# compute normals at the nodes on the surface of solid 3D
function compute_normals!(solid::FEM3D)
    elems   = solid.surface_elems
    xx      = solid.pos
    normals = solid.normals
    #centroids = solid.centroids
    basis_S = solid.basis_S
    #normal  = zeros(3)
    if  typeof(solid.basis_S) <: Quad4 
	N = zeros(4)
	xieta = [0.,0.]
    elseif  typeof(solid.basis_S) <: Tri3
	N= zeros(3)
	xieta = [0.3333333333333,0.3333333333333]
    else
	@error("only supported Quad4 and Tri3 surface elements.")
    end
    d       = Dict{Int, Vector{Vector{Float64}}}()
    # loop over surface elements
    # get the normal for each element
    # store this normal to the nodes of each element using dict 'd'
    @inbounds for ip = 1:size(elems,1)
	elemNodes  =  @view elems[ip,:]  		
	elemNodes0 =        elems[ip,:]  
	coords     =  @view xx[elemNodes]
	detJ       = lagrange_basis!(N, basis_S, xieta, coords)
	if detJ < 0.
	    #  	elemNodes[2] = elemNodes0[3]			  	
	    # elemNodes[3] = elemNodes0[2]
	    # coords       =  @view xx[elemNodes]
	    println("Surface elements with negative Jacobian!!!")
	end
	normal = getNormals(coords,basis_S)
	#normals[ip]   += normal
	# centroids[ip] += 0.333333333333*@SVector[coords[1][1]+coords[2][1]+coords[3][1],
	#                                          coords[1][2]+coords[2][2]+coords[3][2],
	#                                          coords[1][3]+coords[2][3]+coords[3][3],
	#                               ]
	for i in elemNodes
	    if haskey(d,i)
		push!(d[i],copy(normal))
	    else
		d[i] = [copy(normal)]
	    end
	end
    end		
    # now that all normals at all nodes on the surface have been found,
    # comute unique normal = average
    bnd_particles = solid.mesh.node_sets["boundary"]
    #c = 1
    for p in bnd_particles
    	normals_at_p = d[p]
    	normals[p]   = mean(normals_at_p)
        #if p == 41
        #    @printf("Computing normals d[%d]=\t", p)
        #    println(d[p])
        #    @printf("normals[%d] = [%.3e, %.3e, %.3e]\n", p, normals[p][1], normals[p][2], normals[p][3])
        #end
    	#centroids[c] = xx[p]
    	#c += 1
    end
end

export FEM2D, FEM3D, FEMAxis, Rigid2D, move, assign_velocity, fixYNodes, fixNodes, rotate, initializeBasis, compute_normals!


end
