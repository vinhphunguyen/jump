# ----------------------------------------------------------------------
#
#                    ***       JUMP       ***
#                Material Point Method in Julia
#
# Copyright (2020) Vinh Phu Nguyen, phu.nguyen@monash.edu
# Civil Engineering, Monash University
# Clayton VIC 3800, Australia
# This software is distributed under the GNU General Public License.
#
# -----------------------------------------------------------------------

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/jMPM/src")
#
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")


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

#function main()

    # problem parameters
	fGravity      = 0.0
	density       = 1000.0
	youngModulus  = 1000.0
	poissonRatio  = 0.3


    # create the grid of a 1 x 1 square, with 20 x 20 cells
    grid1 =  Grid2D(0.0,0.4,0.0,0.4, 11, 11)
    grid2 =  Grid2D(0.6,1.0, 0.6, 1.0, 11, 11)
    basis = LinearBasis()

    rad     = 0.2
	ppc     = 2
    fOffset = grid1.dx/ppc
    dx      = fOffset
    coords1 = buildParticleForCircle([0.2; 0.2], rad, fOffset)
    coords2 = buildParticleForCircle([0.8; 0.8], rad, fOffset)

    material = ElasticMaterial(youngModulus,poissonRatio,density)

    solid1 = Solid2D(coords1,material)
    solid2 = Solid2D(coords2,material)

    solid1.mass          .*= dx * dx
    solid1.volume        .= dx * dx
    solid1.volumeInitial .= dx * dx
	#
    solid2.mass          .*= dx * dx
    solid2.volume        .= dx * dx
    solid2.volumeInitial .= dx * dx

    v0 = SVector{2,Float64}([0.1  0.1])

    # assign initial velocity for the particles
    assign_velocity(solid1, v0)
    assign_velocity(solid2,-v0)

	buildContactParticlesForCircle(rad, solid1)
	buildContactParticlesForCircle(rad, solid2)

    solids = [solid1,solid2]
	grids  = [grid1,grid2]

    # @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    # @printf("Mass : %+.6e \n", sum(solid1.mass)+sum(solid2.mass))
    # @printf("Vol  : %+.6e \n", sum(solid1.volume)+sum(solid2.volume))
    # @printf("Vol0 : %+.6e \n", sum(solid1.volumeInitial)+sum(solid2.volumeInitial))

    Tf = 3.5 #3.5e-0
    interval=100

    bodyforce = ConstantBodyForce2D(fGravity)

	#output1  = PyPlotOutput(interval,"two-disks-tl/","Two Disks Collision",(4., 4.))
	output2  = OvitoOutput(interval,"two-disks-tl/",["color", "pressure"])
	#fix     =

    problem  = ExplicitDynamicsTL2D(grids,solids,basis,Tf,bodyforce,output2,[])
    algo1    = USL(1e-9)
    #solve(problem,algo1,dtime=0.001)

    plotParticles(problem.output,solids,[grid1.lx, grid1.ly],[grid1.nodeCountX, grid1.nodeCountY],0)
	plotParticles(problem.output,grid1,0)
	plotParticles(problem.output,grid2,1)

    solve(problem,algo1,0.001)

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
	PyPlot.plot(problem.recordTime, c="blue", problem.kinEnergy, "-", label="\$ K \$", linewidth=1.0)
	#PyPlot.hold(true)
	PyPlot.plot(problem.recordTime, c="red", problem.strEnergy, "-", label="\$ U \$", linewidth=1.0)
	PyPlot.plot(problem.recordTime, c="green", problem.kinEnergy + problem.strEnergy, "-", label="\$ K+U \$", linewidth=1.0)
	PyPlot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=8)
	PyPlot.savefig("plot_2Disk_Julia.pdf")

# end
#
# @time main()
