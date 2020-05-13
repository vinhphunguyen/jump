
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

using Solid
using Grid
using Problem
using Output
using Algorithm
using Material
using BodyForce
using Basis
using Fix
using Util

function main()

    # problem parameters
	fGravity      = 0.0
	density       = 1000.0
	youngModulus  = 1000.0
	poissonRatio  = 0.3

    # create the grid of a 1 x 1 square, with 20 x 20 cells
	# and a basis: linear and CPDI-Q4 supported
    grid      =  Grid2D(0,1.5, 0,1.5, 11, 11)
    basis     =  CPDIQ4Basis()

    material = ElasticMaterial(youngModulus,poissonRatio,density,0,0)

    solid1   = Solid2D("disk.msh",material)
    solid2   = Solid2D("disk.msh",material)

    v0 = SVector{2,Float64}([ 0.1  0.1])

    # assign initial velocity for the particles
    assign_velocity(solid1, v0)
    assign_velocity(solid2,-v0)
    # as the mesh was created with the center of the disk at (0,0)
	#move(solid1,SVector{2,Float64}([ 0.2+grid.dx  0.2+grid.dx]))
	#move(solid2,SVector{2,Float64}([ 0.8-grid.dx  0.8-grid.dx]))
	move_cpdi(solid1,SVector{2,Float64}([ 0.45,  0.45]))
	move_cpdi(solid2,SVector{2,Float64}([ 0.45+.6,  .45+.6]))

    solids = [solid1, solid2]

    Tf       = 3. #3.5e-0
    interval = 100
	dtime    = 1e-3

	#output1  = PyPlotOutput(interval,"twodisks-results/","Two Disks Collision",(4., 4.))
	output2  = OvitoOutput(interval,"twodisks-results/",["pressure"])
	fix      = EnergiesFix(solids,"twodisks-results/energies.txt")

    algo1    = USL(1e-9)
    algo2    = MUSL(.99)

	report(grid,solids,dtime)

	#plotParticles(problem.output,solids,[grid.lx, grid.ly],[grid.nodeCountX, grid.nodeCountY],0)
    plotParticles(output2,grid,0)

	#reset_timer!()
    solve_explicit_dynamics_2D(grid,solids,basis,algo2,output2,fix,Tf,dtime)
    #print_timer()
	# plotting energies
    # pyFig_RealTime = PyPlot.figure("MPM 2Disk FinalPlot", figsize=(8/2.54, 4/2.54))
	# PyPlot.clf()
	# pyPlot01 = PyPlot.gca()
	# PyPlot.subplots_adjust(left=0.15, bottom=0.25, right=0.65)
	# pyPlot01[:grid](b=true, which="both", color="gray", linestyle="-", linewidth=0.5)
	# pyPlot01[:set_axisbelow](true)
	# pyPlot01[:set_xlim](0.0, 4.0)
	# pyPlot01[:set_ylim](0.0, 3.0)
	# pyPlot01[:set_xlabel]("time (s)", fontsize=8)
	# pyPlot01[:set_ylabel]("energy (\$\\times 10^{-3}\$ Nm)", fontsize=8)
	# pyPlot01[:set_xticks](collect(0.0:1.0:4.0))
	# pyPlot01[:tick_params](axis="both", which="major", labelsize=8)
	# pyPlot01[:set_yticks](collect(0.0:1.0:3.0))
	# PyPlot.plot(fix.recordTime, c="blue", fix.kinEnergy, "-", label="\$ K \$", linewidth=1.0)
	# #PyPlot.hold(true)
	# PyPlot.plot(fix.recordTime, c="red", fix.strEnergy, "-", label="\$ U \$", linewidth=1.0)
	# PyPlot.plot(fix.recordTime, c="green", fix.kinEnergy + fix.strEnergy, "-", label="\$ K+U \$", linewidth=1.0)
	# PyPlot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=8)
	# #PyPlot.savefig("plot_2Disk_Julia.pdf")

end

@time main()
