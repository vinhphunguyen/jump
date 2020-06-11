
# Phu Nguyen, Monash University
# 20 March, 2020 (Coronavirus outbreak)

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/juMP")
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
using Fix

function main()

    # problem parameters
	fGravity      = 0.0
	density       = 1000.0
	youngModulus  = 1000.0
	poissonRatio  = 0.3


    # create the grid of a 1 x 1 x 1 square, with 10 x 10 x 1 cells
    grid  =  Grid3D(0,1.0, 0,1.0, 0,1.0, 11, 11, 6)
    basis = LinearBasis()
    basis = QuadBsplineBasis()

    fOffset = 0.2/8 # there are 8 material points over the radius (16 MPs)
    # how to calculate fOffset:
    # length of 1 cell = 1/20=0.05 => there are 0.2/0.05=4 cells over the half of the
    # circle. If you want 2 MPs per cell (in x dir.), then there are 8 MPs there.
    ppc = 4
    dx  = fOffset

    coords1 = buildParticleForSphere([0.2; 0.2; 0.5], 0.2, fOffset)
    coords2 = buildParticleForSphere([0.8; 0.8; 0.5], 0.2, fOffset)

    material = ElasticMaterial(youngModulus,poissonRatio,density)

    solid1   = Solid3D(coords1,material)
    solid2   = Solid3D(coords2,material)

    solid1.mass          .*= dx * dx * dx
    solid1.volume        .=  dx * dx * dx
    solid1.volumeInitial .=  dx * dx *dx

    solid2.mass          .*= dx * dx * dx
    solid2.volume        .=  dx * dx * dx
    solid2.volumeInitial .=  dx * dx * dx

    v0 = SVector{3,Float64}([0.1  0.1 0.0])

    # assign initial velocity for the particles
    Solid.assign_velocity(solid1, v0)
    Solid.assign_velocity(solid2,-v0)

    solids = [solid1, solid2]

    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("Mass : %+.6e \n", sum(solid1.mass)+sum(solid2.mass))
    @printf("Vol  : %+.6e \n", sum(solid1.volume)+sum(solid2.volume))
    @printf("Vol0 : %+.6e \n", sum(solid1.volumeInitial)+sum(solid2.volumeInitial))

    Tf      = 3.5 #3.5e-0
    interval= 120
	dtime   = 1e-3

	output2  = OvitoOutput(interval,"results3D/",["pressure", "vx", "vz"])
	fix      = EnergiesFix(solids,"results3D/energies.txt")


    algo2    = MUSL(1.)

    solve_explicit_dynamics_3D(grid,solids,basis,algo2,output2,fix,Tf,dtime)

    pyFig_RealTime = PyPlot.figure("MPM 2Disk FinalPlot", figsize=(8/2.54, 4/2.54))
	PyPlot.clf()
	pyPlot01 = PyPlot.gca()
	PyPlot.subplots_adjust(left=0.15, bottom=0.25, right=0.65)
	pyPlot01[:grid](b=true, which="both", color="gray", linestyle="-", linewidth=0.5)
	pyPlot01[:set_axisbelow](true)
	pyPlot01[:set_xlim](0.0, 5.0)
	pyPlot01[:set_ylim](0.0, 1.0)
	pyPlot01[:set_xlabel]("time (s)", fontsize=8)
	pyPlot01[:set_ylabel]("energy (\$\\times 10^{-3}\$ Nm)", fontsize=8)
	pyPlot01[:set_xticks](collect(0.0:1.0:4.0))
	pyPlot01[:tick_params](axis="both", which="major", labelsize=8)
	pyPlot01[:set_yticks](collect(0.0:0.2:1.0))
	PyPlot.plot(fix.recordTime, c="blue", fix.kinEnergy, "-", label="\$ K \$", linewidth=1.0)
	#PyPlot.hold(true)
	PyPlot.plot(fix.recordTime, c="red", fix.strEnergy, "-", label="\$ U \$", linewidth=1.0)
	PyPlot.plot(fix.recordTime, c="green", fix.kinEnergy + fix.strEnergy, "-", label="\$ K+U \$", linewidth=1.0)
	PyPlot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=8)
	PyPlot.savefig("plot_3Disk_Julia.pdf")

end

@time main()
