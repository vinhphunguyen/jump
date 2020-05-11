
# Phu Nguyen, Monash University
# 20 March, 2020 (Coronavirus outbreak)

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/jMPM/src")
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

function main()

    # problem parameters
	fGravity      = 0.0
	density       = 8940e-12
	youngModulus  = 115e3
	poissonRatio  = 0.31
	yieldStres    = 65.
    alumK         = 10.    # hardening modulus

	length_cyl    = 25.4 #mm
	radius_cyl    = 3.8 #mm


    # create the grid of a 1 x 1 x 1 square, with 10 x 10 x 1 cells
    grid  = Grid3D(14.0, 14.0, 28.0, 20, 20, 40)
    basis = LinearBasis()


    ppc = 4
    fOffset = grid.dx/ppc

    coords1 = buildParticleForCylinder([7.0; 7.0; 0.5*length_cyl], radius_cyl,
	                                   legnth_cyl, fOffset, fOffset, ZAxis)

    material = ElastoPlasticMaterial(youngModulus,poissonRatio,density,yieldStres,alumK,length(coords1))
    solid1   = Solid3D(coords1,material)

    solid1.mass          .*= dx * dx * dx
    solid1.volume        .=  dx * dx * dx
    solid1.volumeInitial .=  dx * dx *dx

    v0 = SVector{3,Float64}([0.0  0. -190000])

    # assign initial velocity for the particles
    assign_velocity(solid1, v0)

    solids = [solid1]

    @printf("Total number of material points: %d \n", solid1.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("Vol  : %+.6e \n", sum(solid1.volume))

    Tf      = 3.5 #3.5e-0
    interval= 120

	fixXForBottom(grid)
	fixYForBottom(grid)
	fixZForBottom(grid)

    bodyforce = ConstantBodyForce3D(fGravity)

	output2  = OvitoOutput(interval,"TaylorBar3D/",["pressure", "vx", "vz"])
	#fix     =

    problem = ExplicitSolidMechanics3D(grid,solids,basis,Tf,bodyforce,output2,[])
    algo1    = USL(1e-9)
    algo2    = MUSL()
    #solve(problem,algo1,dtime=0.001)
    solve(problem,algo1,0.001)

end

@time main()
