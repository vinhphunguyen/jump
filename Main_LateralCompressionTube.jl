# Phu Nguyen, Monash University
# 20 March, 2020 (Coronavirus outbreak)

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/juMP")
# import Gadfly
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")
using BenchmarkTools
using TimerOutputs

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
using Util

function calibrate()
    steelRho  = 7800e-12
    steelE    = 210.0e3
    steelNu   = 0.0
	steelFy   = 210.
	steelK    = 100.
	vel       = 0.1

	grid  = Grid2D(0.,2., 0.,1.,3,2)
	basis = LinearBasis()
	
    # solid creation
    ppc      = 3
	fOffset  = grid.dx/ppc
	coords1  = buildParticleForRectangle([0.5; 0.5], 1.0, grid.ly, fOffset)
	coords2  = buildParticleForRectangle([1.5; 0.5], 1.0, grid.ly, grid.dx/1)


	material1 = ElastoPlasticMaterial(steelE,steelNu,steelRho,steelFy,steelK,length(coords1))
	material2 = RigidMaterial(vel,0.,steelRho)
	

	c_dil     = sqrt((material1.lambda + 2*material1.mu)/material1.density)
	dt        = grid.dx/c_dil
	dtime     = 0.9 * dt

	solid1 = Solid2D(coords1,material1)
	solid2 = Solid2D(coords2,material2)

    solid1.mass          .*=fOffset * fOffset
    solid1.volume        .= fOffset * fOffset
    solid1.volumeInitial .= fOffset * fOffset


    solids = [solid1,solid2]

    # Boundary conditions

    fixXForLeft(grid)
    #fixYForLeft(grid)
    #fixXForRight(grid)

    Tf      = 0.5
    interval= 100

	output2   = OvitoOutput(interval,"tubes/",["pstrain", "vonMises","vy"])
	fix       = StressFix(solids,"tubes/stress-strain.txt")
    algo1     =  MUSL(0.99)
    algo2     =  USL(0.)

	report(grid,solids,dtime)

	#plotParticles_2D(output2,solids,[grid.lx, grid.ly],[grid.nodeCountX, grid.nodeCountY],0)
    plotGrid(output2,grid,0)
    solve_explicit_dynamics_2D(grid,solids,basis,algo2,output2,fix,Tf,dtime)
	

	
end


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
	maxDia = innerDiameter+2*thickness
	grid  = Grid2D(0.,100.0, 0.,maxDia+4, 101, 61)
	basis = LinearBasis()
	basis = QuadBsplineBasis()
    # solid creation
    ppc      = 10
	fOffset  = grid.dx/ppc
	fOffsetR = grid.dx/2
	coords1  = buildParticleForRing([grid.lx/2; grid.ly/2], 0.5*innerDiameter, 0.5*(innerDiameter+2*thickness), fOffset)
	coords2  = buildParticleForRectangle([grid.lx/2; grid.ly-0.5], grid.lx, 1.0, fOffsetR)
	coords3  = buildParticleForRectangle([grid.lx/2; 0.5], grid.lx, 1.0, fOffsetR)

	material1 = ElastoPlasticMaterial(steelE,steelNu,steelRho,steelFy,steelK,length(coords1))
	material2 = RigidMaterial(0.,-v,steelRho)
	material3 = RigidMaterial(0.,0.,steelRho)

	c_dil     = sqrt((material1.lambda + 2*material1.mu)/material1.density)
	dt        = grid.dx/c_dil
	dtime     = 0.2 * dt

	solid1 = Solid2D(coords1,material1)
    solid2 = Solid2D(coords2,material2)
    solid3 = Solid2D(coords3,material3)

    solid1.mass          .*=fOffset * fOffset
    solid1.volume        .= fOffset * fOffset
    solid1.volumeInitial .= fOffset * fOffset


    solids = [solid1, solid2, solid3]

	#bodyforce = ConstantBodyForce2D(fGravity)


    # Boundary conditions

    # fixXForBottom(grid)
    # fixYForBottom(grid)

    Tf      = 1e3
    interval= 1000

	output2   = OvitoOutput(interval,"tubes/",["pstrain", "vonMises","vy"])
	fix       = EmptyFix()
    algo1     =  MUSL(0.99)
    algo2     =  USL(0.)

	report(grid,solids,dtime)

	#plotParticles_2D(output2,solids,[grid.lx, grid.ly],[grid.nodeCountX, grid.nodeCountY],0)
    plotGrid(output2,grid,0)
    solve_explicit_dynamics_2D(grid,solids,basis,algo2,output2,fix,Tf,dtime)

	# reset_timer!()
	# @timeit "1" solve_explicit_dynamics_2D(grid,solids,basis,algo2,output2,fix,Tf,dtime)
	# print_timer()
end

@time calibrate()
