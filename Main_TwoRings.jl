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


# Input file for the two rings collision problem, made of hyperleastic material.
# Solved with the standard MPM method, that is ULMPM with either hat functions, quadratic
# b-splines or cubic b-splines
# Output in folder "rings-ulmpm-results/", with lammps dump files and energies.txt
# Status: negative Jacobian even with cubic splines

push!(LOAD_PATH,"./")
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
using Fix

function main()

    # problem parameters
	fGravity      = 0.0
	density       = 1010e-9
	K             = 121.7e-3         # bulk modulus
	mu            = 26.1e-3          # shear modulus

	E             = 9*K*mu/(3*K+mu)
	nu            = E/2/mu-1.0

	rin           = 30.
	rout          = 40.
	w             = 100.
	l             = 200.

	vel           = 30 # mm/s


    # create the grid  with nx x ny cells
    nx    = 251
    ny    = 251
    grid  =  Grid2D(0, l, 0, w, nx , ny)
    #basis = QuadBsplineBasis()
    basis = CubicBsplineBasis()
    #basis = LinearBasis()

	ppc     = 3
    fOffset = grid.dx/ppc
    dx      = fOffset
    coords1 = buildParticleForRing([l/4; w/2], rin, rout, fOffset)
    coords2 = buildParticleForRing([0.75*l; w/2], rin, rout, fOffset)


    material  = NeoHookeanMaterial(E,nu,density)

	c_dil     = sqrt((material.lambda + 2*material.mu)/material.density)
	dt        = grid.dx/c_dil
	dtime     = 0.6 * dt


    solid1   = Solid2D(coords1,material)
    solid2   = Solid2D(coords2,material)

    solid1.mass          .*= dx * dx
    solid1.volume        .= dx * dx
    solid1.volumeInitial .= dx * dx

    solid2.mass          .*= dx * dx
    solid2.volume        .= dx * dx
    solid2.volumeInitial .= dx * dx

    v0 = SVector{2,Float64}([vel  0.0])

    # assign initial velocity for the particles
    Solid.assign_velocity(solid1, v0)
    Solid.assign_velocity(solid2,-v0)

    solids = [solid1, solid2]

    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("Mass : %+.6e \n", sum(solid1.mass)+sum(solid2.mass))
    @printf("Vol  : %+.6e \n", sum(solid1.volume)+sum(solid2.volume))
    @printf("Vol0 : %+.6e \n", sum(solid1.volumeInitial)+sum(solid2.volumeInitial))
    @printf("Sound vel : %+.6e \n", c_dil)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))

    Tf      = 3.5
    interval= 100


    data               = Dict()
    data["total_time"] = Tf
    data["dt"]         = dtime
    data["time"]       = 0.

    bodyforce = ConstantBodyForce2D([0,0])

	output2  = OvitoOutput(interval,"rings-ulmpm-results/",["vx","sigmaxx"])
	fix      = EnergiesFix(solids,"rings-ulmpm-results/energies.txt")
    algo1    = USL(1e-18)
    algo2    = MUSL(0.99) # 0.99 is the mixing beta PIC-FLIP parameter
    
	plotGrid(output2,grid,0)

	# start solving
    solve_explicit_dynamics_2D(grid,solids,basis,algo1,output2,fix,data)

    # post-processing: plotting energies in time
    pyFig_RealTime = PyPlot.figure("MPM 2Disk FinalPlot", figsize=(8/2.54, 4/2.54))
	PyPlot.clf()
	pyPlot01 = PyPlot.gca()
	PyPlot.subplots_adjust(left=0.15, bottom=0.25, right=0.65)
	pyPlot01[:grid](b=true, which="both", color="gray", linestyle="-", linewidth=0.5)
	pyPlot01[:set_axisbelow](true)
	pyPlot01[:set_xlim](0.0, 100.0)
	#pyPlot01[:set_ylim](0.0, 3.0)
	pyPlot01[:set_xlabel]("time (s)", fontsize=8)
	pyPlot01[:set_ylabel]("energy (\$\\times 10^{-3}\$ Nm)", fontsize=8)
	pyPlot01[:set_xticks](collect(0.0:10.0:100.0))
	pyPlot01[:tick_params](axis="both", which="major", labelsize=8)
	#pyPlot01[:set_yticks](collect(0.0:1.0:3.0))
	PyPlot.plot(problem.recordTime, c="blue", problem.kinEnergy, "-", label="\$ K \$", linewidth=1.0)
	#PyPlot.hold(true)
	PyPlot.plot(problem.recordTime, c="red", problem.strEnergy, "-", label="\$ U \$", linewidth=1.0)
	PyPlot.plot(problem.recordTime, c="green", problem.kinEnergy + problem.strEnergy, "-", label="\$ K+U \$", linewidth=1.0)
	PyPlot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=8)
	PyPlot.savefig("plot_2Disk_Julia.pdf")

end

@time main()
