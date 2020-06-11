
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
	g              = 0
	density        = 7850e-12
	youngModulus1  = 210e3
	youngModulus2  = 210e4
	poissonRatio   = 0.3

    # create the grid of a 1 x 1 square, with 20 x 20 cells
	# and a basis: linear and CPDI-Q4 supported
    grid      =  Grid3D(0,310,0,310,0,20,65,65,3)
    basis     =  LinearBasis()


    material1 = ElasticMaterial(youngModulus1,poissonRatio,density)
    material2 = ElasticMaterial(youngModulus2,poissonRatio,density)

    solid1   = FEM3D("ring1.msh",material1)
    solid2   = FEM3D("ring2.msh",material2)
    
    # as the mesh was created with the center of the disk at (0,0)
	#move(solid1,SVector{2,Float64}([ 0.2+grid.dx  0.2+grid.dx]))
	#move(solid2,SVector{2,Float64}([ 0.8-grid.dx  0.8-grid.dx]))


    alpha = 0.
	Fem.rotate(solid1,alpha)
	Fem.rotate(solid2,alpha)

	Fem.move(solid1,SVector{3,Float64}([155,155,5]))
	Fem.move(solid2,SVector{3,Float64}([155,155,5]))

  
    # Symmetric BCs
	fixYNodes(solid1, "fix")
	fixYNodes(solid2, "fix")

    solids = [solid1 solid2]

    Tf       = 200e-6 #3.5e-0
    interval = 10
	dtime    = 0.4*grid.dx/sqrt(youngModulus2/density)

	#output1  = PyPlotOutput(interval,"twodisks-results/","Two Disks Collision",(4., 4.))
	output2  = VTKOutput(interval,"cylinder-pressure-femp/",["pressure"])
	fix      = DisplacementFemFix(solid1,"cylinder-pressure-femp/",2)
	fix      = EnergiesFix(solids,"cylinder-pressure-femp/energies.txt")

    algo1    = USL(0.)
    algo2    = TLFEM(0.)

    body     = ConstantBodyForce3D(@SVector[0.,-g,0.])

	report(grid,solids,dtime)

    plotGrid(output2,grid,0)
    plotParticles_3D(output2,solids,0)

	#reset_timer!
    solve_explicit_dynamics_femp_3D(grid,solids,basis,body,algo2,output2,fix,Tf,dtime)
    #print_timer()

	# #PyPlot.savefig("plot_2Disk_Julia.pdf")

# end

# @time main()
