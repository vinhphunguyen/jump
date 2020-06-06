
# Phu Nguyen, Monash University
# 20 March, 2020 (Coronavirus outbreak)

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/juMP")
# import Gadfly
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")

# pyFig_RealTime = PyPlot.figure("MPM 2Disk Real-time",
#                                figsize=(8/2.54, 8/2.54), edgecolor="white", facecolor="white")

#include("./Grid.jl")
#include("./Problem.jl")

using Solid
using Grid
using Problem
using Output
using Algorithm
using Material
using BodyForce
using Basis
using Fix

#function main()

    # problem parameters
    density       = 8940e-12
    E             = 115e3
    nu            = 0.31
    A             = 65.
    B             = 356.
    C             = 0.013
    n             = 0.37
    eps0dot       = 1.0

	l0    = 25.4 #mm
	r0    = 3.8 #mm


    # create the grid of a 1 x 1 x 1 square, with 10 x 10 x 1 cells
    grid  = Grid3D(0,14.0,0, 28.0,0,14, 20, 60, 20)
    basis = LinearBasis()


    ppc = 2
    fOffset = grid.dx/ppc

    coords1 = buildParticleForCylinder([7.0;  15; 7.0], r0,
	                                   l0/2, fOffset, fOffset)

    
    material  = JohnsonCookMaterial(E,nu,density,A,B,C,n,eps0dot,grid.dy,length(coords1))
    solid1   = Solid3D(coords1,material)

    solid1.mass          .*= fOffset * fOffset * fOffset
    solid1.volume        .=  fOffset * fOffset * fOffset
    solid1.volumeInitial .=  fOffset * fOffset *fOffset

    v0 = SVector{3,Float64}([0.0 -190000 0])

    # assign initial velocity for the particles
    assign_velocity(solid1, v0)

    solids = [solid1]

    @printf("Total number of material points: %d \n", solid1.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("Vol  : %+.6e \n", sum(solid1.volume))


	#fixXForBottom(grid)
	fixYForBottom(grid)
	#fixZForBottom(grid)


    Tf      = 63e-3
    interval= 500

    c_dil     = sqrt(E/density)
    dt        = grid.dy/c_dil
    dtime     = 0.1 * dt


    output2  = OvitoOutput(interval,"TaylorBar3D/",["pressure", "vy", "pstrain"])
    fix      = EmptyFix()
    algo1    = USL(0.)
    algo2    = MUSL(0.99)

    plotGrid(output2,grid,0)

    plotParticles_3D(output2,solids,[grid.lx, grid.ly, grid.lz],
                        [grid.nodeCountX, grid.nodeCountY, grid.nodeCount],0)

    solve_explicit_dynamics_3D(grid,solids,basis,algo1,output2,fix,Tf,dtime)

# end

# @time main()
