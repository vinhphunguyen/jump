# Phu Nguyen, Monash University
# 20 March, 2020 (Coronavirus outbreak)

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/jMPM/src")
# import Gadfly
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")
using BenchmarkTools

#include("./Grid.jl")
#include("./Problem.jl")

using Material
using Algorithm
using Solid
using Grid
using Problem
using Output
using Fix
using Basis
using BodyForce

function main()
	fGravity  = 0.0

	steelRho  = 7800e-12
    steelE    = 210.0e3
    steelNu   = 0.3
	steelFy   = 310.
	steelK    = 10.

    innerDiameter = 47.9 #mm
	thickness     = 1.48 #mm

	v             = 10e3 #mm/s platen velocity

	# grid creation
	grid  = Grid2D(100.0, 60.0, 51, 21)
	basis = CPDIQ4Basis()
    # solid creation
    ppc      = 10
	fOffset  = grid.dx/ppc
	fOffsetR = grid.dx/2
	coords2  = buildParticleForRectangle([50.0; 60.0-grid.dx/2], 100.0, grid.dx, fOffsetR)
	coords3  = buildParticleForRectangle([50.0; grid.dx/2], 100.0, grid.dx, fOffsetR)

	material1 = ElastoPlasticMaterial(steelE,steelNu,steelRho,steelFy,steelK,7898)
	material2 = RigidMaterial([0.,-v],steelRho)
	material3 = RigidMaterial([0.,0.],steelRho)

	c_dil     = sqrt((material1.lambda + 2*material1.mu)/material1.density)
	dt        = grid.dx/c_dil
	dtime     = 0.2 * dt

	solid1 = Solid2D("tube.msh",material1)
    solid2 = Solid2D(coords2,material2)
    solid3 = Solid2D(coords3,material3)

    solids = [solid1, solid2, solid3]

	bodyforce = ConstantBodyForce2D(fGravity)

	@printf("	Tube,   number of material points: %d \n", solid1.parCount)
	@printf("	Timestep:   %f \n", dtime)

    # Boundary conditions

    # fixXForBottom(grid)
    # fixYForBottom(grid)


    Tf      = 1e3
    interval= 1000

	output2  = OvitoOutput(interval,"tubes-cpdi/",["pstrain", "vonMises"])
	#fix      = DisplacementFix(solids,[29.8 40.],"impact-results/")


	problem   = ExplicitSolidMechanics2D(grid,solids,basis,Tf,bodyforce,output2,[])
    algo1     =  MUSL()
    algo2     =  USL(0.)

	plotParticles(problem.output,solids,[grid.lx, grid.ly],[grid.nodeCountX, grid.nodeCountY],0)
    plotParticles(problem.output,grid,0)
    solve(problem, algo2, dtime)
end

@btime main()
