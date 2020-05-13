# Phu Nguyen, Monash University
# 20 March, 2020 (Coronavirus outbreak)

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/juMP")
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
using Util

#function main()
	steelRho  = 8250e-12
    steelE    = 197.6e3
    steelNu   = 0.3
	steelFy   = 280.
	steelK    = 15.

    innerDiameter = 4.0 #mm
    outerDiameter = 5.0 #mm
	nx            = 2 # number of cells along x-dir
	ny            = 2 # number of cells along y-dir
    l0            = outerDiameter * nx + 10.
    h0            = outerDiameter * ny + 2 # two outer skins of 1mm each
	braze_length  = 1.5

	v             = 10e3 #mm/s platen velocity

	noElemsX      = 81
	noElemsY      = 81

	# grid creation
	grid  = Grid2D(0.,l0,0,h0,noElemsX, noElemsY)
	basis = LinearBasis()
    # solid creation
    ppc      = 6
	fOffset  = grid.dx/ppc
	fOffsetR = grid.dx/2
	coords1  = buildParticleForRingWithBraze([outerDiameter/2; outerDiameter/2], 0.5*innerDiameter, 0.5*outerDiameter, braze_length , fOffset)
    coords   = make_rectangular_pattern(coords1,nx,ny,dx=outerDiameter,dy=outerDiameter)
	coords2  = buildParticleForRectangle([grid.lx/2;0.5],         grid.lx, grid.dx, fOffsetR)
	coords3  = buildParticleForRectangle([grid.lx/2;grid.ly-.5], grid.lx, grid.dx, fOffsetR)


	material1 = ElastoPlasticMaterial(steelE,steelNu,steelRho,steelFy,steelK,length(coords))
	material2 = RigidMaterial(0.,-v,steelRho)
	material3 = RigidMaterial(0.,0.,steelRho)

	c_dil     = sqrt((material1.lambda + 2*material1.mu)/material1.density)
	dt        = grid.dx/c_dil
	dtime     = 0.2 * dt

	solid1    = Solid2D(coords, material1)
    solid2    = Solid2D(coords2,material2)
    solid3    = Solid2D(coords3,material3)

	move(solid1,[5.,1.])

    solids    = [solid1, solid2, solid3]

    Tf       = 1e-8#1e3
    interval = 1000

	output2  = OvitoOutput(interval,"cellular-braze/",["pstrain", "vonMises"])
	fix      = EmptyFix()

    algo1     =  MUSL(.99)
    algo2     =  USL(0.)

	report(grid,solids,dtime)

	plotParticles_2D(output2,solids,[grid.lx, grid.ly],[grid.nodeCountX, grid.nodeCountY],0)
    plotGrid(output2,grid,0)
    #solve_explicit_dynamics_2D(grid,solids,basis,algo2,output2,fix,Tf,dtime)
# end
#
# @time main()
