
# Phu Nguyen, Monash University
# 20 March, 2020 (Coronavirus outbreak)

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/jMPM/src")
# import Gadfly
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")
using BenchmarkTools

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
	fGravity      = 100.0
	density       = 1000.0
	youngModulus  = 1000.0
	poissonRatio  = 0.3


    grid     = Grid2D(1.0, 1.0, 11, 11) # create the grid of a 1 x 1 square, with 20 x 20 cells
    basis    = LinearBasis()
    material = ElasticMaterial(youngModulus,poissonRatio,density)
    solid1   = Solid2D("face.png",1., 1., 0., material) # get particles from black pixels (0)

    solids = [solid1]

	fixXForBottom(grid)
    fixYForBottom(grid)

    @printf("Total number of material points: %d \n", solid1.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)


    Tf = 3.5 #3.5e-0
    interval=100

    bodyforce = ConstantBodyForce2D(fGravity)

	#output1  = PyPlotOutput(interval,"results/","Two Disks Collision",(4., 4.))
	output2  = OvitoOutput(interval,"image-results/",["pressure"])
	#fix     =

    problem = ExplicitSolidMechanics2D(grid,solids,basis,Tf,bodyforce,output2,[])
    algo1    = USL(1e-9)
    algo2    = MUSL()
	algo3    = APIC()
    #solve(problem,algo1,dtime=0.001)
    solve(problem,algo1,0.001)
end

@btime main()
