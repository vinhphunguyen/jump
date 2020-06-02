module Fem

using StaticArrays
using Material
using Mesh
using Grid

###########################################################
# Struct for 2D solid: a set of data for 2D particles
# Createa solid by: solid = Solid2D(coords,mat)
# Createa solid by: solid = Solid2D(coords,mat)
###########################################################
struct FEM2D{T <: MaterialType}
	mass                :: Vector{Float64}
	volume              :: Vector{Float64}

	pos0                :: Vector{SVector{2,Float64}}  # node coords initial
	pos                 :: Vector{SVector{2,Float64}}  # position
	velocity            :: Vector{SVector{2,Float64}}  # velocity
	dU                  :: Vector{SVector{2,Float64}}  # incremental displacements
	fint                :: Vector{SVector{2,Float64}}  # internal forces at FE nodes

	deformationGradient :: Vector{SMatrix{2,2,Float64,4}}  # F, 2x2 matrix
	strain              :: Vector{SMatrix{2,2,Float64,4}}  # stress, 2x2 matrix
	stress              :: Vector{MMatrix{2,2,Float64,4}}  # strain
	

	parCount            :: Int64              # same as element count in the mesh
	nodeCount           :: Int64              # number of nodes in the mesh

	mat                 :: T

	elems               
	# cells               :: Array{MeshCell,1}
	# index               :: Int64

	# particles from a mesh
	function FEM2D(fileName,mat::T) where {T <: MaterialType}
		mesh      = read_GMSH(fileName)
		nodeCount = size(nodes,2)                   # node count
		parCount = size(elems,1)                 # element count
		Identity = SMatrix{2,2}(1, 0, 0, 1)
		F        = fill(Identity,parCount)
		strain   = fill(zeros(2,2),parCount)
		stress   = fill(zeros(2,2),parCount)
		vol      = fill(0.,parCount)
		m        = fill(0.,nodeCount)
		x        = fill(zeros(2),nodeCount)
		nodesX   = fill(zeros(2),nodeCount)
		velo     = fill(zeros(2),nodeCount)
		#println(size(nodes,2))
		#println(parCount)

		for i=1:nodeCount
			nodesX[i] = @SVector [nodes[1,i], nodes[2,i]]
		end


		new{T}(m,vol,nodesX,copy(nodesX),velo,copy(velo),copy(velo),F,
			   strain,stress,parCount,nodeCount,mat,elems)
	end
end

###########################################################
# Struct for 3D solid: a set of data for 2D particles
# Createa solid by: solid = Solid2D(coords,mat)
# Createa solid by: solid = Solid2D(coords,mat)
###########################################################
struct FEM3D{T <: MaterialType}
	mass                :: Vector{Float64}
	volume              :: Vector{Float64}

	pos0                :: Vector{SVector{3,Float64}}  # node coords initial
	pos                 :: Vector{SVector{3,Float64}}  # position
	velocity            :: Vector{SVector{3,Float64}}  # velocity
	dU                  :: Vector{SVector{3,Float64}}  # incremental displacements
	fint                :: Vector{SVector{3,Float64}}  # internal forces at FE nodes
	fbody               :: Vector{SVector{3,Float64}}  # external forces  due to gravity at FE nodes

	deformationGradient :: Vector{SMatrix{3,3,Float64,9}}  # F, 2x2 matrix
	strain              :: Vector{SMatrix{3,3,Float64,9}}  # stress, 2x2 matrix
	stress              :: Vector{MMatrix{3,3,Float64,9}}  # strain
	

	parCount            :: Int64              # same as element count in the mesh
	nodeCount           :: Int64              # number of nodes in the mesh

	mat                 :: T

	elems               :: Array{Int64,2}
	mesh                :: FEMesh
	fixedNodesX         :: Vector{Int64}
	fixedNodesY         :: Vector{Int64}
	fixedNodesZ         :: Vector{Int64}

	# particles from a mesh
	function FEM3D(fileName,mat::T) where {T <: MaterialType}
		mesh      = read_GMSH(fileName)
		nodeCount = length(mesh.nodes)                       # node count
		parCount  = length(mesh.element_sets["All"])         # element count
		Identity  = SMatrix{3,3}(1, 0, 0,0, 1, 0,0, 0, 1)
		F         = fill(Identity,parCount)
		strain    = fill(zeros(3,3),parCount)
		stress    = fill(zeros(3,3),parCount)
		vol       = fill(0.,parCount)
		m         = fill(0.,nodeCount)
		x         = fill(zeros(3),nodeCount)
		nodesX    = fill(zeros(3),nodeCount)
		elems     = Array{Array{Int64,1},1}(undef,0)   # element nodes
		velo      = fill(zeros(3),nodeCount)
		#println(size(nodes,2))
		#println(parCount)

        # convert from mesh.nodes to our traditional data structure
		for i=1:nodeCount			
			nodesX[i] = mesh.nodes[i]
		end

        # convert from mesh.elements to our traditional data structure
        volumetricElems = collect(mesh.element_sets["All"])
		for i=1:parCount			
			push!(elems,mesh.elements[volumetricElems[i]])
		end

		fixX       = fill(0,nodeCount)
		fixY       = fill(0,nodeCount)
		fixZ       = fill(0,nodeCount)

		new{T}(m,vol,nodesX,copy(nodesX),velo,copy(velo),copy(velo),copy(velo),F,
			   strain,stress,parCount,nodeCount,mat,vcat(map(x->x', elems)...), mesh, fixX,fixY,fixZ)
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

export FEM2D, FEM3D, move, assign_velocity, fixYNodes


end