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

module Solid
    using Statistics
    using LinearAlgebra  # if not yet installed, in REPL, do import Pkg and Pkd.add("LinearAlgebra")
    using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")
	using Images
	using Printf


    using Material
	using Mesh
	using Grid

	@enum Direction begin
       XAxis
       YAxis
	   ZAxis
    end

	struct Solid1D
		mass                :: Vector{Float64}
		volumeInitial       :: Vector{Float64}
		volume              :: Vector{Float64}

		pos                 :: Vector{Float64}  # position
		pos0                :: Vector{Float64}  # old position
		velocity            :: Vector{Float64}  # velocity

		deformationGradient :: Vector{Float64}  # F, x2 matrix

		strain              :: Vector{Float64}  # stress,
		stress              :: Vector{Float64}  # strain

		parCount            :: Int64   # cannot be UInt because size(A) and loop ares Int

		mat                 :: MaterialType

		function Solid1D(coords::Vector{Float64},mat::MaterialType)
			parCount = length(coords)
			F        = fill(1.,parCount)
			strain   = fill(1.,parCount)
			stress   = fill(1.,parCount)
			m        = fill(mat.density,parCount)
			vol      = fill(0.,parCount)
			vol0     = fill(0.,parCount)
			velo     = fill(0.,parCount)
			new(m,vol,vol0,coords,copy(coords),velo,F,strain,stress,parCount,mat)
		end
	end

    ###########################################################
    # Struct for 2D solid: a set of data for 2D particles
	# Createa solid by: solid = Solid2D(coords,mat)
	# Createa solid by: solid = Solid2D(coords,mat)
    ###########################################################
	struct Solid2D{T <: MaterialType}
		mass                :: Vector{Float64}
		volumeInitial       :: Vector{Float64}
		volume              :: Vector{Float64}

		temp                :: Vector{Float64}             # temperature

		pos0                :: Vector{SVector{2,Float64}}  # position
		pos                 :: Vector{SVector{2,Float64}}  # position
		velocity            :: Vector{SVector{2,Float64}}  # velocity

        deformationGradient :: Vector{SMatrix{2,2,Float64,4}}  # F, 2x2 matrix
		strain              :: Vector{SMatrix{2,2,Float64,4}}  # stress, 2x2 matrix
		stress              :: Vector{SMatrix{2,2,Float64,4}}  # strain
		gradVelo            :: Vector{MMatrix{2,2,Float64,4}}  # velocity gradient
		Cmat                :: Vector{MMatrix{2,2,Float64,4}}  # APIC BpDp

		parCount            :: Int64                       # number of total particle in the solid

		mat                 :: T

        # for CPDIs
		nodes               :: Vector{SVector{2,Float64}}  # position
		elems

		rigid               :: Bool                          # rigid body or

		color               :: Vector{Float64} # used to classify surface particles
		damage              :: Vector{Float64} # for contact in TLMPM
		gradDamage          :: Vector{SVector{2,Float64}}
		crackForce          :: Vector{Float64} # for contact in TLMPM
		contactPar          :: Vector{Int64}   # used to classify surface particles
		radius              :: Vector{Float64} # for contact in TLMPM

        #Point(x::T, y::T) where {T<:Real} = Point{T}(x,y);

        # if particle coords are given
        function Solid2D(coords::Vector{SVector{2,Float64}},mat::T) where {T <: MaterialType}
        	parCount = length(coords)
			Identity = SMatrix{2,2}(1, 0, 0, 1)
			F        = fill(Identity,parCount)
			strain   = fill(zeros(2,2),parCount)
			stress   = fill(zeros(2,2),parCount)
			gradVel  = fill(zeros(2,2),parCount)
			Cmat     = fill(zeros(2,2),parCount)
			m        = fill(mat.density,parCount)
			vol      = zeros(parCount)
			vol0     = zeros(parCount)
			col      = zeros(parCount)
			dam      = zeros(parCount)
			#x        = fill(zeros(2),parCount)
			velo     = fill(zeros(2),parCount)
			rigid    = false
			if ( typeof(mat) <: RigidMaterial )
				rigid = true
			end
			return new{T}(m,vol,vol0,copy(vol0), copy(coords),coords,velo,F,strain,stress,gradVel,
			    Cmat,parCount,mat,fill(zeros(2),1),0, rigid,col, dam, copy(velo), copy(dam))
        end

        #if particles are CPDI, from a mesh
		function Solid2D(nodes,elems,mat::T) where {T <: MaterialType}
			parCount = size(elems,1)
			Identity = SMatrix{2,2}(1, 0, 0, 1)
			F        = fill(Identity,parCount)
			strain   = fill(zeros(2,2),parCount)
			stress   = fill(zeros(2,2),parCount)
			gradVel  = fill(zeros(2,2),parCount)
			Cmat     = fill(zeros(2,2),parCount)
			m        = fill(mat.density,parCount)
			vol      = zeros(parCount)
			vol0     = zeros(parCount)
			x        = fill(zeros(2),parCount)
			nodesX   = fill(zeros(2),size(nodes,2))
			velo     = fill(zeros(2),parCount)

			for e = 1:parCount
			    coord =  nodes[:,elems[e,:]]
				a     = 0.5*( coord[1,1]*coord[2,2]  - coord[1,2]*coord[2,1]
						    + coord[1,2]*coord[2,3]  - coord[1,3]*coord[2,2]
						    + coord[1,3]*coord[2,4]  - coord[1,4]*coord[2,3]
						    + coord[1,4]*coord[2,1]  - coord[1,1]*coord[2,4])
			    vol[e]  = a;
			    vol0[e] = a;
			    m[e]    = a*mat.density;
			    x[e]    = vec(mean(coord,dims=2)); # center of each element=particle
			end

			for i=1:size(nodes,2)
				nodeX[i] = @SVector [nodes[1,i], nodes[2,i]]
			end

			new{T}(m,vol,vol0,copy(x),x,velo,F,strain,stress,gradVel,
			    Cmat,parCount,mat,nodesX,elems)
		end
        # particles from a mesh: only Q4 for the moment
		function Solid2D(fileName,mat::T) where {T <: MaterialType}
			nodes,elems = loadGMSH(fileName)


			parCount = size(elems,1)
			Identity = SMatrix{2,2}(1, 0, 0, 1)
			F        = fill(Identity,parCount)
			strain   = fill(zeros(2,2),parCount)
			stress   = fill(zeros(2,2),parCount)
			gradVel  = fill(zeros(2,2),parCount)
			Cmat     = fill(zeros(2,2),parCount)
			m        = fill(mat.density,parCount)
			vol      = zeros(parCount)
			vol0     = zeros(parCount)
			x        = fill(zeros(2),parCount)
			nodesX   = fill(zeros(2),size(nodes,2))
			velo     = fill(zeros(2),parCount)
			#println(size(nodes,2))
			#println(parCount)
			
			for i=1:size(nodes,2)
				nodesX[i] = @SVector [nodes[1,i], nodes[2,i]]
			end

			rigid    = false
			if ( typeof(mat) <: RigidMaterial )
				rigid = true
			end


			# check negative Jacobian issue
			N         = zeros(4)#@SVector [0,0,0,0]
			@inbounds for ip = 1:parCount
				elemNodes  =  @view elems[ip,:]  
				elemNodes0 =        elems[ip,:]  
				coords     =  @view nodesX[elemNodes]
		    
		
			    detJ  = lagrange_basis!(N, Quad4(), [0., 0.], coords)
			    if detJ < 0.
			  	  elemNodes[2] = elemNodes0[4]			  	
			  	  elemNodes[4] = elemNodes0[2]
			  	  println(elemNodes)
			    end

			    coord = nodes[1:2,elems[ip,:]] # from gmsh, nodes are 3D
			    a     = 0.5*( coord[1,1]*coord[2,2]  - coord[1,2]*coord[2,1]
				            + coord[1,2]*coord[2,3]  - coord[1,3]*coord[2,2]
							+ coord[1,3]*coord[2,4]  - coord[1,4]*coord[2,3]
							+ coord[1,4]*coord[2,1]  - coord[1,1]*coord[2,4])
			    vol[ip]  = a
			    vol0[ip] = a
			    m[ip]    = a*mat.density
			    x[ip]    = vec(mean(coord,dims=2)) # center of each element=particle
			end


			new{T}(m,vol,vol0,copy(x),x,velo,F,strain,stress,gradVel,Cmat,parCount,mat,
			    nodesX,elems,rigid)
		end

		# particles from an image 'fileName'
		# lx = physical length of the image in x-dir
		# all pixels are used!!!
		function Solid2D(fileName,lx, ly, phaseId, material::T) where {T <: MaterialType}
			img      = load(fileName)               # read the image
			imgg     = Gray.(img)                   # convert to gray scale
			mat      = convert(Array{Float64}, imgg) # convert to float64 matrix
			dx       = lx/size(mat,2)
			dy       = ly/size(mat,1)

            xVec     = Vector{SVector{2,Float64}}(undef,0)  # position
			for j = 1:size(mat,2)
			   for i = 1:size(mat,1)
				    val = mat[i,j]
					if abs(val - phaseId) < 1e-10
						x = (j-1)*dx + dx/2;
                        y = (i-1)*dy + dy/2;
						push!(xVec,[x,y])
					end
			   end
		    end
			parCount = length(xVec)
			Identity = SMatrix{2,2}(1, 0, 0, 1)
			F        = fill(Identity,parCount)
			strain   = fill(zeros(2,2),parCount)
			stress   = fill(zeros(2,2),parCount)
			gradVel  = fill(zeros(2,2),parCount)
			Cmat     = fill(zeros(2,2),parCount)
			m        = fill(dx*dy*material.density,parCount)
			vol      = fill(dx*dy,parCount)
			vol0     = fill(dx*dy,parCount)
			velo     = fill(zeros(2),parCount)

			new{T}(m,vol,vol0,copy(xVec),xVec,velo,F,strain,stress,gradVel,
			    Cmat,parCount,material,0,0)
		end
   end

   #Solid2D(coords::Vector{SVector{2,Float64}},mat::T) where {T <: MaterialType} = Solid2D{T}(coords,mat)

   ##################################################
   ###      Solid3D
   ##################################################

   struct Solid3D{T<:MaterialType}
	   mass                :: Vector{Float64}
	   volumeInitial       :: Vector{Float64}
	   volume              :: Vector{Float64}

	   pos                 :: Vector{SVector{3,Float64}}  # position
	   velocity            :: Vector{SVector{3,Float64}}  # velocity

	   deformationGradient :: Vector{SMatrix{3,3,Float64,9}}  # F, 3x3 matrix
	   strain              :: Vector{SMatrix{3,3,Float64,9}}  # strain, 3x3 matrix
	   stress              :: Vector{SMatrix{3,3,Float64,9}}  # stress

	   parCount            :: Int64

	   mat                 :: T

	   function Solid3D(coords::Vector{SVector{3,Float64}},mat::T) where {T <: MaterialType}
		   parCount = length(coords)
		   Identity = SMatrix{3,3}(1, 0, 0,  0, 1, 0, 0, 0, 1)
		   F        = fill(Identity,parCount)
		   dF       = fill(Identity,parCount)
		   strain   = fill(zeros(3,3),parCount)
		   stress   = fill(zeros(3,3),parCount)
		   m        = fill(mat.density,parCount)
		   vol      = fill(0,parCount)
		   vol0     = fill(0,parCount)
		   #x        = fill(zeros(2),parCount)
		   velo     = fill(zeros(3),parCount)
		   return new{T}(m,vol0,vol,coords,velo,F,strain,stress,parCount,mat)
	   end
  end

	function toXArray(solid::Solid2D)
		return [(solid.pos[i])[1] for i = 1:solid.parCount]
	end

	function toYArray(solid::Solid2D)
		return [(solid.pos[i])[2] for i = 1:solid.parCount]
	end

	function buildParticleForSegment(length, fOffset)
		coord   = Vector{Float64}(undef,0)

		for fx = 0.5*fOffset:fOffset:length-0.5*fOffset
			push!(coord, fx)
		end

		return(coord)
	end

   """
       buildParticleForCircle(fCenter, fRadius, fOffset)
   generate particles for simple geometries
   """

   function buildParticleForCircle(fCenter::Array{Float64}, fRadius::Float64, fOffset::Float64)
		thisMaterialDomain = Vector{SVector{2,Float64}}(undef,0)
		coord              = SVector{2,Float64}

		fRadius = floor(fRadius/fOffset) * fOffset	#just in case radius is not a multiple of offset
		for fy in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
			for fx in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
				if(fx^2 + fy^2 < fRadius^2)
					coord = [fCenter[1] + fx; fCenter[2] + fy]
					push!(thisMaterialDomain, coord)
				end
			end
		end

		return(thisMaterialDomain)
	end

	function buildParticleForRing(fCenter::Array{Float64}, innerRad, outerRad, fOffset::Float64)
	   thisMaterialDomain = Vector{SVector{2,Float64}}(undef,0)
	   coord              = SVector{2,Float64}

	   fRadius = floor(outerRad/fOffset) * fOffset	#just in case radius is not a multiple of offset
	   for fy in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
		   for fx in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
			   if (fx^2 + fy^2 < outerRad^2) && (fx^2 + fy^2 > innerRad^2)
				   coord = [fCenter[1] + fx; fCenter[2] + fy]
				   push!(thisMaterialDomain, coord)
			   end
		   end
	   end

	   return(thisMaterialDomain)
   end

	function buildParticleForSphere(fCenter::Array{Float64}, fRadius::Float64, fOffset::Float64)
	   thisMaterialDomain = Vector{SVector{3,Float64}}(undef,0)
	   coord              = SVector{3,Float64}

	   fRadius = floor(fRadius/fOffset) * fOffset	#just in case radius is not a multiple of offset
	   for fz in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
		   for fy in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
			   for fx in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
				   if(fx^2 + fy^2 +fz^2 < fRadius^2)
					   coord = [fCenter[1] + fx; fCenter[2] + fy; fCenter[3]+ fz]
					   push!(thisMaterialDomain, coord)
				   end
			   end
		   end
       end

	   return(thisMaterialDomain)
   end

   # ring with braze of length l0
   # cellular materials
   function buildParticleForRingWithBraze(fCenter::Array{Float64},
	                                      innerRad, outerRad, l0, fOffset, fOffsetB)
	  thisMaterialDomain = Vector{SVector{2,Float64}}(undef,0)
	  coord              = SVector{2,Float64}

	  fRadius = floor(outerRad/fOffset) * fOffset	#just in case radius is not a multiple of offset
	  for fy in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
		  for fx in -fRadius+0.5*fOffset:fOffset:+fRadius-0.5*fOffset
			  if (fx^2 + fy^2 < outerRad^2) && (fx^2 + fy^2 > innerRad^2)
				  coord = [fCenter[1] + fx; fCenter[2] + fy]
				  push!(thisMaterialDomain, coord)
			  end
		  end
	  end

	  # handling the vertical braze, which is a rectagle
	  AB      = 0.5*l0
	  R       = outerRad
	  x       = sqrt(R^2-AB^2)
	  CB      = R - x
      count = 0
	  # println(AB)
	  # println(CB)
	  for fy in -AB+0.5*fOffsetB:fOffsetB:+AB-0.5*fOffsetB
		  for fx in -CB+0.5*fOffsetB:fOffsetB:+CB-0.5*fOffsetB
			  # compute the coords in global coords
			  coord = [fCenter[1] + R + fx; fCenter[2] + fy]
			  # coords with origin at the ring center
			  dx1   = coord[1] - fCenter[1]
			  dx2   = coord[2] - fCenter[2]

			  if (dx1^2 + dx2^2 > R^2)
				  # coords with origin at 0' (center of the right ring)
	  			  xx = coord[1] - fCenter[1] - 2*R
	  			  yy = coord[2] - fCenter[2]
				  if ( xx^2 + yy^2 > R^2)
				    push!(thisMaterialDomain, coord)
					count += 1
			      end
			  end
		  end
	  end
	  # for the horizontal braze_length
	  for fy in -CB+0.5*fOffsetB:fOffsetB:+CB-0.5*fOffsetB
		  for fx in -AB+0.5*fOffsetB:fOffsetB:+AB-0.5*fOffsetB
			  # compute the coords in global coords
			  coord = [fCenter[1] + fx; fCenter[2] + R + fy]
			  # coords with origin at the ring center
			  dx1   = coord[1] - fCenter[1]
			  dx2   = coord[2] - fCenter[2]

			  if (dx1^2 + dx2^2 > R^2)
				  # coords with origin at 0' (center of the right ring)
	  			  xx = coord[1] - fCenter[1]
	  			  yy = coord[2] - fCenter[2] - 2*R
				  if ( xx^2 + yy^2 > R^2)
				    push!(thisMaterialDomain, coord)
					count += 1
			      end
			  end
		  end
	  end
	  @printf("Number of braze material points: %d \n", count)

	  return(thisMaterialDomain)
	 end

   ###############################################################
   # buildParticleForCylinder, y-axis
   ###############################################################
   function buildParticleForCylinder(fCenter::Array{Float64}, fRadius, fLength,
	                                 fOffsetX, fOffsetZ)
	 thisMaterialDomain = Vector{SVector{3,Float64}}(undef,0)
	 coord              = SVector{3,Float64}

	 fRadius = floor(fRadius/fOffsetX) * fOffsetX	#just in case radius is not a multiple of offset
	 for fy in -fLength+0.5*fOffsetZ:fOffsetZ:+fLength-0.5*fOffsetZ
		 for fz in -fRadius+0.5*fOffsetX:fOffsetX:+fRadius-0.5*fOffsetX
			 for fx in -fRadius+0.5*fOffsetX:fOffsetX:+fRadius-0.5*fOffsetX
				 if(fx^2 + fz^2 < fRadius^2)
					 coord = [fCenter[1] + fx; fCenter[2] + fy; fCenter[3] + fz]
					 push!(thisMaterialDomain, coord)
				 end
			 end
		 end
     end

	 return(thisMaterialDomain)
 end

	"""
        buildParticleForRectangle(fCenter, fRadius, fOffset)
    generate particles for simple geometries
    """
	function buildParticleForRectangle(fCenter::Vector{Float64}, fWidth::Float64, fHeight::Float64, fOffset::Float64)
		thisMaterialDomain = Vector{SVector{2,Float64}}(undef,0)
		coord              = SVector{2,Float64}

        fWidth	= floor(fWidth/fOffset) * fOffset	#just in case width is not a multiple of offset
		fHeight	= floor(fHeight/fOffset) * fOffset	#just in case height is not a multiple of offset

		for fy = -0.5*fHeight+0.5*fOffset:fOffset:+0.5*fHeight-0.5*fOffset
			for fx = -0.5*fWidth+0.5*fOffset:fOffset:+0.5*fWidth-0.5*fOffset
				coord = [fCenter[1] + fx; fCenter[2] + fy]
				push!(thisMaterialDomain, coord)
			end
		end

		return(thisMaterialDomain)
	end

		"""
        buildParticleForBlock(fCenter, fRadius, fOffset)
    generate particles for simple geometries
    """
	function buildParticleForBlock(fCenter::Vector{Float64}, fWidth, fHeight, fThickness, fOffset)
		thisMaterialDomain = Vector{SVector{3,Float64}}(undef,0)
		coord              = SVector{3,Float64}

        fWidth	= floor(fWidth/fOffset) * fOffset	#just in case width is not a multiple of offset
		fHeight	= floor(fHeight/fOffset) * fOffset	#just in case height is not a multiple of offset

		for fz = -0.5*fThickness+0.5*fOffset:fOffset:+0.5*fThickness-0.5*fOffset
			for fy = -0.5*fHeight+0.5*fOffset:fOffset:+0.5*fHeight-0.5*fOffset
				for fx = -0.5*fWidth+0.5*fOffset:fOffset:+0.5*fWidth-0.5*fOffset
					coord = [fCenter[1] + fx; fCenter[2] + fy; fCenter[3] + fz]
					push!(thisMaterialDomain, coord)
				end
			end
		end

		return(thisMaterialDomain)
	end

	function buildParticleForRectangleWithANotch(fCenter::Vector{Float64},
		fWidth::Float64, fHeight::Float64, fOffset::Float64, x1, x2, ymax)
		thisMaterialDomain = Vector{SVector{2,Float64}}(undef,0)
		coord              = SVector{2,Float64}

		fWidth	= floor(fWidth/fOffset) * fOffset	#just in case width is not a multiple of offset
		fHeight	= floor(fHeight/fOffset) * fOffset	#just in case height is not a multiple of offset

		for fy = -0.5*fHeight+0.5*fOffset:fOffset:+0.5*fHeight-0.5*fOffset
			for fx = -0.5*fWidth+0.5*fOffset:fOffset:+0.5*fWidth-0.5*fOffset
				coord = [fCenter[1] + fx; fCenter[2] + fy]
				if ( coord[1] > x1 && coord[1] < x2 && coord[2] < ymax ) continue end
				push!(thisMaterialDomain, coord)
			end
		end

		return(thisMaterialDomain)
	end


    # assign an initial velocity to all particles of solid
	function assign_velocity(solid,v0)
       velo = solid.velocity
       @inbounds for p = 1 : solid.parCount
          velo[p] = v0
       end
	end

	function move_cpdi(solid,dx)
	   x1 = solid.nodes
	   @inbounds for p = 1 : length(x1)
		  x1[p] += dx
	   end
	   x2 = solid.pos
	   @inbounds for p = 1 : solid.parCount
		  x2[p] += dx
	   end
	end

	function move(solid,dx)
       x = solid.pos
       @inbounds for p = 1 : solid.parCount
          x[p] += dx
       end
	end

	function rotate(solid,alpha)
	   xx = solid.pos
	   beta = deg2rad(alpha)
	   @inbounds for p = 1 : solid.parCount
		  xp = xx[p][1]
		  yp = xx[p][2]
		  xq = xp*cos(beta)-yp*sin(beta)
		  yq = yp*cos(beta)+xp*sin(beta)
		  xx[p] = setindex(xx[p],xq,1)
		  xx[p] = setindex(xx[p],yq,2)
	   end
	end

    # from a given coords (which is actually a geometry)
	# make a rectanglular pattern of nx x ny objects with
	# distance in x dir = dx, and in y-dir = dy
	# coded for cellular structure modeling
	function make_rectangular_pattern(coords,nx,ny;dx=0.,dy=0.)
		res      = Vector{SVector{2,Float64}}(undef,0)
		ptsCount = length(coords)
		for j=1:ny
			for i=1:nx
              for p=1:ptsCount
				  # coords of the original point
				  x0 = coords[p][1]
				  y0 = coords[p][2]
				  # coors of the generated point
				  x  = x0 + (i-1) * dx
				  y  = y0 + (j-1) * dy
				  push!(res,@SVector[x,y])
			  end
			end
		end
		return res
	end

	# function doCellParticleInteraction(solid::Solid2D,grid::Grid2D)
   #    xx     = solid.pos
   #    pElems = repeat(0:0, inner=solid.parCount)
   #    for p=1:solid.parCount
   #        x = xx[p][1] - grid.xmin
   #        y = xx[p][2] - grid.ymin
   #        e = floor(Int64,x/grid.dx) + 1 + (grid.nodeCountX-1)*floor(Int64,y/grid.dy)
   #        pElems[p] = e
   #    end
   #
   #     mpoints = Vector{Vector{Int64}}(undef,0)
   #
   #     for e=1:(grid.nodeCountX-1)*(grid.nodeCountY-1)
   #        id  = findall(pElems->pElems==e,vec(pElems))
   #        push!(mpoints,id)
   #     end
   #     return mpoints
   # end
   #
   # function buildContactParticlesForCircle(fRadius::Float64, solid)
	#    XX = solid.pos0
	#    @inbounds for ip = 1:solid.parCount
  	# 	 xp = XX[ip][1]
  	# 	 yp = XX[ip][2]
	# 	 if abs(xp^2 + yp^2 - fRadius^2) < 1e-10
	# 		 push!(solid.contactPar,ip)
	# 	 end
	#  end
   # end

	export buildParticleForCircle, buildParticleForRing, buildParticleForRectangle, buildParticleForRectangleWithANotch, buildParticleForSegment, rotate, buildParticleForBlock,
	buildParticleForSphere, buildParticleForRingWithBraze, buildParticleForCylinder, toXArray,toYArray,assign_velocity, move, move_cpdi, make_rectangular_pattern
	export Solid1D, Solid2D, Solid3D
end
