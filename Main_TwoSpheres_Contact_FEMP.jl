
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
	fGravity      = 0.0
	density       = 1000.0
	youngModulus  = 1000.0
	poissonRatio  = 0.3

    # create the grid of a 1 x 1 square, with 20 x 20 cells
	# and a basis: linear and CPDI-Q4 supported
    grid      =  Grid3D(0,1.,0,1.,0,1., 31, 31,31)
    basis     =  LinearBasis()


    solid1   = FEM3D("sphere-1.msh")
    solid2   = FEM3D("sphere-1.msh")


    material = ElasticMaterial(youngModulus,poissonRatio,density,solid1.parCount)
    #material = NeoHookeanMaterial(youngModulus,poissonRatio,density)

    v0 = SVector{3,Float64}([ 0.1  0.1 0.0])

    # assign initial velocity for the particles
    Fem.assign_velocity(solid1, v0)
    Fem.assign_velocity(solid2,-v0)
    # as the mesh was created with the center of the disk at (0,0)
	#move(solid1,SVector{2,Float64}([ 0.2+grid.dx  0.2+grid.dx]))
	#move(solid2,SVector{2,Float64}([ 0.8-grid.dx  0.8-grid.dx]))
	Fem.move(solid1,SVector{3,Float64}([ 0.2,   0.2,0.5]))
	Fem.move(solid2,SVector{3,Float64}([ 0.2+0.6,0.8,0.5]))

    solids = [solid1, solid2]
    mats   = [material,material]

    Tf       = 3.5 #3.5e-0
    interval = 100
	dtime    = 1e-3

	#output1  = PyPlotOutput(interval,"twodisks-results/","Two Disks Collision",(4., 4.))
	output2  = VTKOutput(interval,"two-spheres-femp/",["pressure"])
	fix      = EnergiesFix(solids,"two-spheres-femp/energies.txt")

    algo2    = USL(1e-9)
    algo1    = TLFEM(1e-9,1.)
    body     = ConstantBodyForce3D([0.,0.,0.])

	report(grid,solids,dtime)

    plotGrid(output2,grid,0)
    plotParticles_3D(output2,solids,mats,0)

	#reset_timer!
    fric = 0.

    data                    = Dict()
    data["total_time"]      = Tf
    data["time"]            = 0.
    data["dt"]              = dtime
    data["friction"]        = fric

    solve_explicit_dynamics_femp_3D_Contact(grid,solids,mats,basis,body,algo1,output2,fix,data)
    #print_timer()
	# plotting energies
    pyFig_RealTime = PyPlot.figure("MPM 2Disk FinalPlot", figsize=(8/2.54, 4/2.54))
	PyPlot.clf()
	pyPlot01 = PyPlot.gca()
	PyPlot.subplots_adjust(left=0.15, bottom=0.25, right=0.65)
	pyPlot01[:grid](b=true, which="both", color="gray", linestyle="-", linewidth=0.5)
	pyPlot01[:set_axisbelow](true)
	pyPlot01[:set_xlim](0.0, 4.0)
	pyPlot01[:set_ylim](0.0, 3.0)
	pyPlot01[:set_xlabel]("time (s)", fontsize=8)
	pyPlot01[:set_ylabel]("energy (\$\\times 10^{-3}\$ Nm)", fontsize=8)
	pyPlot01[:set_xticks](collect(0.0:1.0:4.0))
	pyPlot01[:tick_params](axis="both", which="major", labelsize=8)
	pyPlot01[:set_yticks](collect(0.0:1.0:3.0))
	PyPlot.plot(fix.recordTime, c="blue", fix.kinEnergy, "-", label="\$ K \$", linewidth=1.0)
	#PyPlot.hold(true)
	PyPlot.plot(fix.recordTime, c="red", fix.strEnergy, "-", label="\$ U \$", linewidth=1.0)
	PyPlot.plot(fix.recordTime, c="green", fix.kinEnergy + fix.strEnergy, "-", label="\$ K+U \$", linewidth=1.0)
	PyPlot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=8)
	# #PyPlot.savefig("plot_2Disk_Julia.pdf")

# end

# @time main()
