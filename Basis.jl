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

# Module Basis:
# for all supported basis functions
module Basis

using  Grid
using  Solid
using  Printf
using  StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")

include("Bsplines.jl")

export BasisType,LinearBasis, QuadBsplineBasis, CPDIQ4Basis
export getShapeAndGradient, getShapeAndGradientTL, getShapeFunctions, getShapeFuncs, initialise
export getAdjacentGridPoints, quad_bspline_type1, quad_bspline_type2,
       quad_bspline_type3, quad_bspline_type4

abstract type BasisType end

struct LinearBasis      <: BasisType      end
struct CPDIQ4Basis      <: BasisType      end
struct QuadBsplineBasis <: BasisType      end

#################################################################
# getAdjacentGridPoints: 1D
# get indices of points adjacent to point 'xp'
# indices : one dimensional indices
function getAdjacentGridPoints(nearPts,xp,grid::Grid1D, basis::LinearBasis)
	iBottomLeft_i  = floor(Int64,xp * grid.dxI + 1.)
	nearPts       .= [iBottomLeft_i, iBottomLeft_i + 1]
	#println(thisAdjacentGridPoints)
end
#################################################################
# getAdjacentGridPoints: 2D
# get indices of points adjacent to point 'xp'
# nearPts : one dimensional indices
function getAdjacentGridPoints(nearPts,xp,grid::Grid2D,basis::LinearBasis)
	fLength_Cell_x = grid.dxI
	fLength_Cell_y = grid.dyI

	iBottomLeft_i  = floor(Int64,(xp[1]-grid.xmin) * fLength_Cell_x + 1.)
	iBottomLeft_j  = floor(Int64,(xp[2]-grid.ymin) * fLength_Cell_y + 1.)

	if(iBottomLeft_j < 1 || iBottomLeft_j > grid.nodeCountY)
		@printf("Index out of bounds: j: %d \n", iBottomLeft_j)
		@printf("xp[2]: %e \n", xp[2])
	end

	iIndex   = index2DTo1D(iBottomLeft_i, iBottomLeft_j,   grid.nodeCountX, grid.nodeCountY)
	nearPts  .= @SVector [iIndex, iIndex+1, iIndex+grid.nodeCountX,iIndex+1+grid.nodeCountX ]
end

# used in CPDI basis functions evaluation
function getAdjacentGridPoints(xp,grid::Grid2D,basis::LinearBasis)
	fLength_Cell_x = grid.dxI
	fLength_Cell_y = grid.dyI

	iBottomLeft_i  = floor(Int64,(xp[1]-grid.xmin) * fLength_Cell_x + 1.)
	iBottomLeft_j  = floor(Int64,(xp[2]-grid.ymin) * fLength_Cell_y + 1.)

	if(iBottomLeft_j < 1 || iBottomLeft_j > grid.nodeCountY)
		@printf("Index out of bounds: j: %d \n", iBottomLeft_j)
		@printf("xp[2]: %e \n", xp[2])
	end

	iIndex   = index2DTo1D(iBottomLeft_i, iBottomLeft_j,   grid.nodeCountX, grid.nodeCountY)
	return [iIndex, iIndex+1, iIndex+grid.nodeCountX,iIndex+1+grid.nodeCountX ]
end

#################################################################
# getAdjacentGridPoints: 3D
# get indices of points adjacent to point 'xp'
# indices : one dimensional indices
function getAdjacentGridPoints(nearPts,xp,grid::Grid3D, basis::LinearBasis)
	fLength_Cell_x = grid.dxI
	fLength_Cell_y = grid.dyI
	fLength_Cell_z = grid.dzI

	iBottomLeft_i  = floor(Int64,xp[1] * fLength_Cell_x + 1.)
	iBottomLeft_j  = floor(Int64,xp[2] * fLength_Cell_y + 1.)
	iBottomLeft_k  = floor(Int64,xp[3] * fLength_Cell_z + 1.)

	if(iBottomLeft_j < 1 || iBottomLeft_j > grid.nodeCountY)
		@printf("Index out of bounds: j: %d \n", iBottomLeft_j)
		@printf("xp[2]: %e \n", xp[2])
	end

	iIndex = index3DTo1D(iBottomLeft_i, iBottomLeft_j, iBottomLeft_k,
	                     grid.nodeCountX, grid.nodeCountY, grid.nodeCountZ)


	nearPts .= @SVector [iIndex, iIndex+1, iIndex+grid.nodeCountX, iIndex+grid.nodeCountX+1,
	                     iIndex+grid.nodeCountXY, iIndex+1+grid.nodeCountXY, iIndex+grid.nodeCountX+grid.nodeCountXY, iIndex+grid.nodeCountX+1+grid.nodeCountXY]
end
#################################################################
# getAdjacentGridPoints: 1D, quadratic bsplines
# get indices of points adjacent to point 'xp'
# indices : one dimensional indices
function getAdjacentGridPoints(nearPts,xp,grid::Grid1D, basis::QuadBsplineBasis)
	iBottomLeft_i  = floor(Int64,xp * grid.dxI + 1.)
	num            = xp * grid.dxI
	whole          = floor(num)
	rem            = num - whole
	#println(iBottomLeft_i)
	if     iBottomLeft_i == 1
	  nearPts .= @SVector[iBottomLeft_i, iBottomLeft_i + 1, iBottomLeft_i + 2]
    elseif iBottomLeft_i == grid.nodeCount-1
	  nearPts .= @SVector[iBottomLeft_i - 1, iBottomLeft_i, iBottomLeft_i + 1]
	else
	  if rem < 0.5
	     nearPts .= @SVector[iBottomLeft_i-1, iBottomLeft_i, iBottomLeft_i+1]
      else
	     nearPts .= @SVector[iBottomLeft_i, iBottomLeft_i+1, iBottomLeft_i+2]
      end
    end
	#println(nearPts)
end

#################################################################
# getShapeAndGradient: 1D, linear basis
# outputs: nearPoints, funcs and and ders
function getShapeAndGradient(nearPoints::Vector{Int64}, funcs::Vector{Float64},
						 ders::Vector{Float64},p::Int64, grid::Grid1D, solid,basis::LinearBasis)
	xp = solid.pos[p]
	getAdjacentGridPoints(nearPoints,xp,grid,basis)

	dxI = grid.dxI

	@inbounds for i = 1:2
		index      = nearPoints[i]
		v2Distance = xp - grid.pos[index]

		Nx   = 1.0 - abs(v2Distance) * dxI
		#if (Nx < 0.0) Nx = 0.0 end, no need
		dNdx = -sign(v2Distance) * dxI

		funcs[i]  = Nx
		ders[i]   = dNdx
	end
	return 2
end

#################################################################
# getShapeAndGradient: 2D, linear basis
# outputs: nearPoints, funcs and and ders
function getShapeAndGradient(nearPoints::Vector{Int64}, funcs::Vector{Float64},
                             ders::Matrix{Float64},p::Int64, grid::Grid2D, solid,basis::LinearBasis)
	xp          = solid.pos[p][1]
	yp          = solid.pos[p][2]
	getAdjacentGridPoints(nearPoints,solid.pos[p],grid,basis)

	dxI = grid.dxI
	dyI = grid.dyI

	@inbounds for i = 1:4
		index      = nearPoints[i]
		xip        = xp - grid.pos[index][1]
		yip        = yp - grid.pos[index][2]

		Nx          = 1.0 - abs(xip) * dxI
		Ny          = 1.0 - abs(yip) * dyI

		funcs[i]    = Nx * Ny
		ders[1,i]   = -Ny*sign(xip) * dxI
		ders[2,i]   = -Nx*sign(yip) * dyI
	end
	return 4
end

function getShapeAndGradientTL(nearPoints::Vector{Int64}, funcs::Vector{Float64},
                               ders::Matrix{Float64},p::Int64, grid::Grid2D, solid,basis::LinearBasis)
	xp          = solid.pos0[p]
	getAdjacentGridPoints(nearPoints,xp,grid,basis)

	dxI = grid.dxI
	dyI = grid.dyI

	@inbounds for i = 1:4
		index      = nearPoints[i]
		xip        = xp[1] - grid.pos[index][1]
		yip        = xp[2] - grid.pos[index][2]

		Nx = 1.0 - abs(xip) * dxI
		Ny = 1.0 - abs(yip) * dyI

		# if (Nx < 0.0) Nx = 0.0 end
		# if (Ny < 0.0) Ny = 0.0 end

		dNdx = -Ny*sign(xip) * dxI
		dNdy = -Nx*sign(yip) * dyI

		funcs[i]    = Nx * Ny
		ders[1,i]   =   dNdx
		ders[2,i]   = 	dNdy
	end

	return 4
end

#################################################################
# getShapeAndGradient: 3D, linear basis
# outputs: nearPoints, funcs and and ders
function getShapeAndGradient(nearPoints::Vector{Int64}, funcs::Vector{Float64},
                             ders::Matrix{Float64},p::Int64, grid::Grid3D, solid,basis::LinearBasis)
	xp = solid.pos[p]
	getAdjacentGridPoints(nearPoints,xp,grid,basis)

	dxI = grid.dxI
	dyI = grid.dyI
	dzI = grid.dzI

	@inbounds for i = 1:8
		index      = nearPoints[i]
		v2Distance = xp - grid.pos[index]

		Nx = 1.0 - abs(v2Distance[1]) * dxI
		Ny = 1.0 - abs(v2Distance[2]) * dyI
		Nz = 1.0 - abs(v2Distance[3]) * dzI

		dNdx = -Ny*Nz*sign(v2Distance[1]) * dxI
		dNdy = -Nx*Nz*sign(v2Distance[2]) * dyI
		dNdz = -Nx*Ny*sign(v2Distance[3]) * dzI

		funcs[i]    = Nx * Ny * Nz
		ders[1,i]   = dNdx
		ders[2,i]   = dNdy
		ders[3,i]   = dNdz
	end
	return 8
end

#################################################################
# getShapeAndGradient: 1D, modified quadratic bspline basis
# outputs: nearPoints, funcs and and ders
function getShapeAndGradient(nearPoints::Vector{Int64}, funcs::Vector{Float64},
						 ders::Vector{Float64},p::Int64, grid::Grid1D, solid, basis::QuadBsplineBasis						 )
	xp          = solid.pos[p]
	getAdjacentGridPoints(nearPoints,xp,grid,basis)

	dxI = grid.dxI

	@inbounds for i = 1:3
		index      = nearPoints[i]
		r          = (xp - grid.pos[index]) * dxI

		if     index == 1 | index == grid.nodeCount
			N,dN = quad_bspline_type1(r)
		elseif index == 2
			N,dN = quad_bspline_type2(r)
		elseif index == grid.nodeCount - 1
			N,dN = quad_bspline_type4(r)
		else
			N,dN = quad_bspline_type3(r)
		end
	    funcs[i] = N
		ders[i]  = dN * dxI
	end
	println(sum(funcs))
	println(sum(ders))
	return 3
end

#################################################################
# getShapeAndGradient: 2D, modified quadratic bspline basis
# outputs: nearPoints, funcs and and ders
function getShapeAndGradient(nearPoints::Vector{Int64}, funcs::Vector{Float64},
						 ders::Matrix{Float64},p::Int64, grid::Grid2D, solid,basis::QuadBsplineBasis)
	xp          = solid.pos[p][1]
	yp          = solid.pos[p][2]

    # inverse of grid spacing
	dxI = grid.dxI
	dyI = grid.dyI


	numx            = (xp - grid.xmin) * grid.dxI
	wholex          = floor(numx)
	remx            = numx - wholex

	numy            = (yp - grid.ymin) * grid.dyI
	wholey          = floor(numy)
	remy            = numy - wholey


	iBottomLeft_i  = floor(Int64,numx + 1.)
	iBottomLeft_j  = floor(Int64,numy + 1.)

	#println(iBottomLeft_i)
	if     iBottomLeft_i == 1
	  nearPointsX = @SVector [iBottomLeft_i, iBottomLeft_i + 1, iBottomLeft_i + 2]
    elseif iBottomLeft_i == grid.nodeCountX-1
	  nearPointsX = @SVector [iBottomLeft_i - 1, iBottomLeft_i, iBottomLeft_i + 1]
	else
	  if remx < 0.5
	     nearPointsX = @SVector [iBottomLeft_i-1, iBottomLeft_i, iBottomLeft_i+1]
      else
	     nearPointsX = @SVector [iBottomLeft_i, iBottomLeft_i+1, iBottomLeft_i+2]
      end
    end

	if     iBottomLeft_j == 1
	  nearPointsY = @SVector [iBottomLeft_j, iBottomLeft_j + 1, iBottomLeft_j + 2]
    elseif iBottomLeft_j == grid.nodeCountY-1
	  nearPointsY = @SVector [iBottomLeft_j - 1, iBottomLeft_j, iBottomLeft_j + 1]
	else
	  if remy < 0.5
		 nearPointsY = @SVector [iBottomLeft_j-1, iBottomLeft_j, iBottomLeft_j+1]
	  else
		 nearPointsY = @SVector [iBottomLeft_j, iBottomLeft_j+1, iBottomLeft_j+2]
	  end
	end

    index = 1
	for i = 1:3
		id     = nearPointsX[i]
		r      = (xp - grid.pos[id][1]) * dxI
        if     id == 1 | id == grid.nodeCountX
			Nx,dNx = quad_bspline_type1(r)
		elseif id == 2
			Nx,dNx = quad_bspline_type2(r)
		elseif id == grid.nodeCountX - 1
			Nx,dNx = quad_bspline_type4(r)
		else
			Nx,dNx = quad_bspline_type3(r)
		end
		for j = 1:3
			jd     = nearPointsY[j]
			r      = (yp - (jd-1)*grid.dy) * dyI
			if     jd == 1 | jd == grid.nodeCountY
				Ny,dNy = quad_bspline_type1(r)
			elseif jd == 2
				Ny,dNy = quad_bspline_type2(r)
			elseif jd == grid.nodeCountY - 1
				Ny,dNy = quad_bspline_type4(r)
			else
				Ny,dNy = quad_bspline_type3(r)
			end
			nearPoints[index] = index2DTo1D(id,jd,grid.nodeCountX, grid.nodeCountY)
		    funcs[index]      = Nx  * Ny
			ders[1,index]     = dNx * dxI * Ny
			ders[2,index]     = dNy * dyI * Nx
			index += 1
		end
    end
	# println(sum(funcs))
	# println(sum(ders))
	return 9
end


#################################################################
# getShapeAndGradient: 3D, modified quadratic bspline basis
# outputs: nearPoints, funcs and and ders
function getShapeAndGradient(nearPoints::Vector{Int64}, funcs::Vector{Float64},
						     ders::Matrix{Float64},p::Int64, grid::Grid3D, solid,basis::QuadBsplineBasis)
	xp          = solid.pos[p][1]
	yp          = solid.pos[p][2]
	zp          = solid.pos[p][3]

    # inverse of grid spacing
	dxI = grid.dxI
	dyI = grid.dyI
	dzI = grid.dzI


	numx            = xp * grid.dxI
	wholex          = floor(numx)
	remx            = numx - wholex

	numy            = yp * grid.dyI
	wholey          = floor(numy)
	remy            = numy - wholey

	numz            = zp * grid.dzI
	wholez          = floor(numz)
	remz            = numz - wholez

	iBottomLeft_i  = floor(Int64,numx + 1.)
	iBottomLeft_j  = floor(Int64,numy + 1.)
	iBottomLeft_k  = floor(Int64,numz + 1.)

	#println(iBottomLeft_i)
	if     iBottomLeft_i == 1
	  nearPointsX = @SVector[iBottomLeft_i, iBottomLeft_i + 1, iBottomLeft_i + 2]
    elseif iBottomLeft_i == grid.nodeCountX-1
	  nearPointsX = @SVector[iBottomLeft_i - 1, iBottomLeft_i, iBottomLeft_i + 1]
	else
	  if remx < 0.5
	     nearPointsX = @SVector[iBottomLeft_i-1, iBottomLeft_i, iBottomLeft_i+1]
      else
	     nearPointsX = @SVector[iBottomLeft_i, iBottomLeft_i+1, iBottomLeft_i+2]
      end
    end

	if     iBottomLeft_j == 1
	  nearPointsY = @SVector[iBottomLeft_j, iBottomLeft_j + 1, iBottomLeft_j + 2]
    elseif iBottomLeft_j == grid.nodeCountY-1
	  nearPointsY = @SVector[iBottomLeft_j - 1, iBottomLeft_j, iBottomLeft_j + 1]
	else
	  if remy < 0.5
		 nearPointsY = @SVector[iBottomLeft_j-1, iBottomLeft_j, iBottomLeft_j+1]
	  else
		 nearPointsY = @SVector[iBottomLeft_j, iBottomLeft_j+1, iBottomLeft_j+2]
	  end
	end

	if     iBottomLeft_k == 1
	  nearPointsZ = @SVector[iBottomLeft_k, iBottomLeft_k + 1, iBottomLeft_k + 2]
    elseif iBottomLeft_j == grid.nodeCountY-1
	  nearPointsZ = @SVector[iBottomLeft_k - 1, iBottomLeft_k, iBottomLeft_k + 1]
	else
	  if remz < 0.5
		 nearPointsZ = @SVector[iBottomLeft_k-1, iBottomLeft_k, iBottomLeft_k+1]
	  else
		 nearPointsZ = @SVector[iBottomLeft_k, iBottomLeft_k+1, iBottomLeft_k+2]
	  end
	end

    index = 1
	for i = 1:3
		id     = nearPointsX[i]
		r      = (xp - grid.pos[id][1]) * dxI
        if     id == 1 | id == grid.nodeCountX
			Nx,dNx = quad_bspline_type1(r)
		elseif id == 2
			Nx,dNx = quad_bspline_type2(r)
		elseif id == grid.nodeCountX - 1
			Nx,dNx = quad_bspline_type4(r)
		else
			Nx,dNx = quad_bspline_type3(r)
		end
		for j = 1:3
			jd     = nearPointsY[j]
			r      = (yp - (jd-1)*grid.dy) * dyI
			if     jd == 1 | jd == grid.nodeCountY
				Ny,dNy = quad_bspline_type1(r)
			elseif jd == 2
				Ny,dNy = quad_bspline_type2(r)
			elseif jd == grid.nodeCountY - 1
				Ny,dNy = quad_bspline_type4(r)
			else
				Ny,dNy = quad_bspline_type3(r)
			end
			for k = 1:3
				kd     = nearPointsZ[j]
				r      = (zp - (kd-1)*grid.dz) * dzI
				if     kd == 1 | kd == grid.nodeCountZ
					Nz,dNz = quad_bspline_type1(r)
				elseif kd == 2
					Nz,dNz = quad_bspline_type2(r)
				elseif kd == grid.nodeCountZ - 1
					Nz,dNz = quad_bspline_type4(r)
				else
					Nz,dNz = quad_bspline_type3(r)
				end
			end
			nearPoints[index] = index3DTo1D(id,jd,kd,grid.nodeCountX, grid.nodeCountY,grid.nodeCountZ)
		    funcs[index]      = Nx  * Ny * Nz
			ders[1,index]     = dNx * dxI * Ny * Nz
			ders[2,index]     = dNy * dyI * Nx * Nz
			ders[3,index]     = dNz * dzI * Nx * Ny
			index += 1
		end
    end
	# println(sum(funcs))
	# println(sum(ders))
	return 27
end


#################################################################
# getShapeFunctions: 1D,  LinearBasis
# outputs: nearPoints, funcs and and ders
function getShapeFunctions(nearPoints::Vector{Int64}, funcs::Vector{Float64},
                       p::Int64, grid::Grid1D, solid,basis::LinearBasis)
	xp = solid.pos[p]
	getAdjacentGridPoints(nearPoints,xp,grid,basis)

	dxI = grid.dxI

	@inbounds for i = 1:2
		index      = nearPoints[i]
		v2Distance = xp - grid.pos[index]

		Nx = 1.0 - abs(v2Distance) * dxI

		if (Nx < 0.0) Nx = 0.0 end
		funcs[i]    = Nx
	end
	return 2
end


#################################################################
# getShapeFunctions: 2D,  LinearBasis
# outputs: nearPoints, funcs and and ders
function getShapeFunctions(nearPoints::Vector{Int64}, funcs::Vector{Float64},
                           p::Int64, grid::Grid2D, solid,basis::LinearBasis)
	xp = solid.pos[p]
	getAdjacentGridPoints(nearPoints,xp,grid,basis)

	dxI = grid.dxI
	dyI = grid.dyI

	@inbounds for i = 1:4
		index      = nearPoints[i]
		v2Distance = xp - grid.pos[index]

		Nx = 1.0 - abs(v2Distance[1]) * dxI
		Ny = 1.0 - abs(v2Distance[2]) * dyI

		funcs[i]    = Nx * Ny
	end
	return 4
end

#################################################################
# getShapeAndGradient: 2D, modified quadratic bspline basis
# outputs: nearPoints, funcs and and ders
function getShapeFunctions(nearPoints::Vector{Int64}, funcs::Vector{Float64},
						   p::Int64, grid::Grid2D, solid,basis::QuadBsplineBasis)
	xp          = solid.pos[p][1]
	yp          = solid.pos[p][2]

    # inverse of grid spacing
	dxI = grid.dxI
	dyI = grid.dyI


	numx            = (xp - grid.xmin) * grid.dxI
	wholex          = floor(numx)
	remx            = numx - wholex

	numy            = (yp - grid.ymin) * grid.dyI
	wholey          = floor(numy)
	remy            = numy - wholey


	iBottomLeft_i  = floor(Int64,numx + 1.)
	iBottomLeft_j  = floor(Int64,numy + 1.)

	#println(iBottomLeft_i)
	if     iBottomLeft_i == 1
	  nearPointsX = @SVector[iBottomLeft_i, iBottomLeft_i + 1, iBottomLeft_i + 2]
    elseif iBottomLeft_i == grid.nodeCountX-1
	  nearPointsX = @SVector[iBottomLeft_i - 1, iBottomLeft_i, iBottomLeft_i + 1]
	else
	  if remx < 0.5
	     nearPointsX = @SVector[iBottomLeft_i-1, iBottomLeft_i, iBottomLeft_i+1]
      else
	     nearPointsX = @SVector[iBottomLeft_i, iBottomLeft_i+1, iBottomLeft_i+2]
      end
    end

	if     iBottomLeft_j == 1
	  nearPointsY = @SVector[iBottomLeft_j, iBottomLeft_j + 1, iBottomLeft_j + 2]
    elseif iBottomLeft_j == grid.nodeCountY-1
	  nearPointsY = @SVector[iBottomLeft_j - 1, iBottomLeft_j, iBottomLeft_j + 1]
	else
	  if remy < 0.5
		 nearPointsY = @SVector[iBottomLeft_j-1, iBottomLeft_j, iBottomLeft_j+1]
	  else
		 nearPointsY = @SVector[iBottomLeft_j, iBottomLeft_j+1, iBottomLeft_j+2]
	  end
	end

    index = 1
	for i = 1:3
		id     = nearPointsX[i]
		r      = (xp - grid.pos[id][1]) * dxI
        if     id == 1 | id == grid.nodeCountX
			Nx,dNx = quad_bspline_type1(r)
		elseif id == 2
			Nx,dNx = quad_bspline_type2(r)
		elseif id == grid.nodeCountX - 1
			Nx,dNx = quad_bspline_type4(r)
		else
			Nx,dNx = quad_bspline_type3(r)
		end
		for j = 1:3
			jd     = nearPointsY[j]
			r      = (yp - (jd-1)*grid.dy) * dyI
			if     jd == 1 | jd == grid.nodeCountY
				Ny,dNy = quad_bspline_type1(r)
			elseif jd == 2
				Ny,dNy = quad_bspline_type2(r)
			elseif jd == grid.nodeCountY - 1
				Ny,dNy = quad_bspline_type4(r)
			else
				Ny,dNy = quad_bspline_type3(r)
			end
			nearPoints[index] = index2DTo1D(id,jd,grid.nodeCountX, grid.nodeCountY)
		    funcs[index]      = Nx  * Ny
			index += 1
		end
    end
	# println(sum(funcs))
	# println(sum(ders))
	return 9
end

# compute linear basis and grads for all nodes I at particle 'p'
# ders: 2x4 matrix, each col is for each node
function getShapeFuncs(nearPoints::Vector{Int64}, funcs::Vector{Float64},xp, grid::Grid2D, solid,basis::LinearBasis)
	# earPoints = getAdjacentGridPoints(xp,grid)
	# funcs      = zeros(4)
	# ders       = zeros(2,4)

	getAdjacentGridPoints(nearPoints,xp,grid,basis)

	dxI = grid.dxI
	dyI = grid.dyI

	@inbounds for i = 1:4
		index      = nearPoints[i]
		v2Distance = xp - grid.pos[index]

		Nx = 1.0 - abs(v2Distance[1]) * dxI
		Ny = 1.0 - abs(v2Distance[2]) * dyI

		if (Nx < 0.0) Nx = 0.0 end
		if (Ny < 0.0) Ny = 0.0 end

		funcs[i]    = Nx * Ny
	end
end

function getShapeIp(xp, xI,grid::Grid2D, basis::LinearBasis)
	v2Distance = xp - xI

	Nx = 1.0 - abs(v2Distance[1]) * grid.dxI
	Ny = 1.0 - abs(v2Distance[2]) * grid.dyI

	if (Nx < 0.0) Nx = 0.0 end
	if (Ny < 0.0) Ny = 0.0 end

	return Nx * Ny
end

#################################################################
# getShapeAndGradient: 2D, CPDI-Q4
# outputs: nearPoints, funcs and and ders

function getShapeAndGradient(nearPoints::Vector{Int64}, funcs::Vector{Float64},
                         ders::Matrix{Float64},p::Int64, grid::Grid2D, solid,basis::CPDIQ4Basis)
	 nodeIds = solid.elems[p,:]
	 corners = solid.nodes[:,nodeIds]

	 x11     = corners[1,1]
	 x12     = corners[2,1]
	 x21     = corners[1,2]
	 x22     = corners[2,2]
	 x31     = corners[1,3]
	 x32     = corners[2,3]
	 x41     = corners[1,4]
	 x42     = corners[2,4]

	 Vp     = 0.5*   ( x11*x22  - x21*x12 +
	                   x21*x32  - x31*x22 +
					   x31*x42  - x41*x32 +
					   x41*x12  - x11*x42 )

	 # function and gradient weights
	 c1   = (x21-x11)*(x42-x12) - (x22-x12)*(x41-x11)
	 c2   = (x21-x11)*(x32-x22) - (x22-x12)*(x31-x21)
	 c3   = (x31-x41)*(x42-x12) - (x32-x42)*(x41-x11)
	 c4   = (x31-x41)*(x32-x22) - (x32-x42)*(x31-x21)

	 wf   = (1/(36*Vp))*[4*c1+2*c2+2*c3+c4 2*c1+4*c2+c3+2*c4
	                       c1+2*c2+2*c3+4*c4 2*c1+c2+4*c3+2*c4]

	 wg   = zeros(2,4)

	 wg[:,1] = [x22-x42 x41-x21]
	 wg[:,2] = [x32-x12 x11-x31]
	 wg[:,3] = -wg[:,1]
	 wg[:,4] = -wg[:,2]

	 wg      *= 1/(2*Vp)

	 # nodes I where phi_I(xp) are non-zero
	 nearPts = Vector{Int64}(undef,16)
	 @inbounds for c=1:4
	     nearPts[(c-1)*4+1:c*4] = getAdjacentGridPoints(corners[:,c],grid)
	 end
	 nearPts   = unique(nearPts)
	 x         = findall(iszero, nearPts)
	 if length(x) > 0 @error("Zero index: %d \n", p) end
	 nodeCount = length(nearPts)
	 nearPoints[1:nodeCount] = nearPts

	 # compute phi_I(xp) and first derivatives

	 funcs .= 0.
	 ders  .= 0.

	 @inbounds for i=1:nodeCount
	     xI = grid.pos[nearPoints[i]]
	     @inbounds for c=1:4
	         N           = getShapeIp(corners[:,c], xI,grid)
	         funcs[i]   += wf[c]  *N
	         ders[:,i]  += wg[:,c]*N
	     end
	 end
	 return nodeCount
end

#################################################################
# getShapeFunctions: 2D, CPDI-Q4
# outputs: nearPoints, funcs and and ders
function getShapeFunctions(nearPoints::Vector{Int64}, funcs::Vector{Float64},
                         p::Int64, grid::Grid2D, solid,basis::CPDIQ4Basis)
	 nodeIds = solid.elems[p,:]
	 corners = solid.nodes[:,nodeIds]

	 x11     = corners[1,1]
	 x12     = corners[2,1]
	 x21     = corners[1,2]
	 x22     = corners[2,2]
	 x31     = corners[1,3]
	 x32     = corners[2,3]
	 x41     = corners[1,4]
	 x42     = corners[2,4]

	 Vp     = 0.5*   ( x11*x22  - x21*x12 +
	                   x21*x32  - x31*x22 +
					   x31*x42  - x41*x32 +
					   x41*x12  - x11*x42 )

	 # function and gradient weights
	 c1   = (x21-x11)*(x42-x12) - (x22-x12)*(x41-x11)
	 c2   = (x21-x11)*(x32-x22) - (x22-x12)*(x31-x21)
	 c3   = (x31-x41)*(x42-x12) - (x32-x42)*(x41-x11)
	 c4   = (x31-x41)*(x32-x22) - (x32-x42)*(x31-x21)

	 wf   = (1/(36*Vp))*[4*c1+2*c2+2*c3+c4 2*c1+4*c2+c3+2*c4
	                       c1+2*c2+2*c3+4*c4 2*c1+c2+4*c3+2*c4]
	# nodes I where phi_I(xp) are non-zero
		nearPts = Vector{Int64}(undef,16)
		@inbounds for c=1:4
			nearPts[(c-1)*4+1:c*4] = getAdjacentGridPoints(corners[:,c],grid)
		end
		nearPts   = unique(nearPts)
		nodeCount = length(nearPts)
		nearPoints[1:nodeCount] = nearPts

	 # compute phi_I(xp) and first derivatives

	 funcs .= 0.

	 @inbounds for i=1:nodeCount
	     xI = grid.pos[nearPoints[i]]
	     @inbounds for c=1:4
	         N           = getShapeIp(corners[:,c], xI,grid)
	         funcs[i]   += wf[c]  *N
	     end
	 end

	 return nodeCount
end

function initialise(grid::Grid1D,basis::LinearBasis)
	nearPoints = [0,0]
	funcs      = [0.,0.]
	ders       = [0.,0.]
	return (nearPoints, funcs, ders)
end

function initialise(grid::Grid2D,basis::LinearBasis)
	nearPoints    = repeat(0:0, inner=4) #[0, 0, 0, 0]
	funcs         = [0., 0., 0., 0.]
	ders          = [0. 0. 0. 0.; 0. 0. 0. 0.]
	return nearPoints, funcs, ders
end

function initialise(grid::Grid3D, basis::LinearBasis)
	nearPoints    = [0, 0, 0, 0, 0, 0, 0, 0]
	funcs         = zeros(8)
	ders          = zeros(3,8)
	return (nearPoints, funcs, ders)
end

function initialise(grid::Grid1D,basis::QuadBsplineBasis)
	nearPoints = [0,0,0]
	funcs      = zeros(3)
	ders       = zeros(3)
	return (nearPoints, funcs, ders)
end

function initialise(grid::Grid2D,basis::QuadBsplineBasis)
	nearPoints = repeat(0:0, inner=9)
	funcs      = zeros(9)
	ders       = zeros(2,9)
	return (nearPoints, funcs, ders)
end

function initialise(grid::Grid3D,basis::QuadBsplineBasis)
	nearPoints = repeat(0:0, inner=27)
	funcs      = zeros(27)
	ders       = zeros(3,27)
	return (nearPoints, funcs, ders)
end

function initialise(grid::Grid2D,basis::CPDIQ4Basis)
	nearPoints = repeat(0:0, inner=16)
	funcs      = zeros(16)
	ders       = zeros(2,16)
	return (nearPoints, funcs, ders)
end



end
