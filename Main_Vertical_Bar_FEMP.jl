
# Phu Nguyen, Monash University
# 20 March, 2020 (Coronavirus outbreak)

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/juMP")
# import Gadfly
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")
using TimerOutputs
# pyFig_RealTime = PyPlot.figure("MPM 2Disk Real-time",
#                                figsize=(8/2.54, 8/2.54), edgecolor="white", facecolor="white")

#include("./Grid.jl")
#include("./Problem.jl")

using Fem
using Solid
using Grid
using FemMPM
using Output
using Algorithm
using Material
using BodyForce
using Basis
using Fix
using Util

#function main()

    # problem parameters
	fGravity      = 1e6
	density       = 1050e-12
	youngModulus  = 1.0
	poissonRatio  = 0.3

    # create the grid of a 1 x 1 square, with 20 x 20 cells
	# and a basis: linear and CPDI-Q4 supported
    grid      =  Grid3D(0,2000,0,3500,0,2000,5,8,5)
    basis     =  LinearBasis()


    material = NeoHookeanMaterial(youngModulus,poissonRatio,density)

    solid1   = FEM3D("bar.msh",material)
    
    # as the mesh was created with the center of the disk at (0,0)
	#move(solid1,SVector{2,Float64}([ 0.2+grid.dx  0.2+grid.dx]))
	#move(solid2,SVector{2,Float64}([ 0.8-grid.dx  0.8-grid.dx]))
	Fem.move(solid1,SVector{3,Float64}([2000/4,2490,2000/4]))

    #Fem.assign_velocity(solid1, SVector{3,Float64}([0. -1000. 0.0 ]))

	fixYForTop(grid)
	fixYForBottom(grid)

    solids = [solid1]

    Tf       = 0.25 #3.5e-0
    interval = 2
	dtime    = 0.1*grid.dx/sqrt(youngModulus/density)

	#output1  = PyPlotOutput(interval,"twodisks-results/","Two Disks Collision",(4., 4.))
	output2  = VTKOutput(interval,"vertical-bar-femp/",["pressure"])
	fix      = DisplacementFemFix(solid1,"vertical-bar-femp/",2)

    algo1    = USL(0.)
    algo1    = TLFEM(0.)

	report(grid,solids,dtime)

    #plotGrid(output2,grid)
    #plotParticles_3D(output2,solids,0)

	#reset_timer!
    solve_explicit_dynamics_femp_3D(grid,solids,basis,algo1,output2,fix,Tf,dtime)
    #print_timer()

	# #PyPlot.savefig("plot_2Disk_Julia.pdf")

# end

# @time main()
