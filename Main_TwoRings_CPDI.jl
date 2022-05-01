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

# Input file for the two rubber rings collision problem
# Solved with the CPDI-Q4 method, 
# Output in folder "cantilever-box-femp/", with lammps dump files and energies.txt

push!(LOAD_PATH,"./")
#
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")


#include("./Grid.jl")
#include("./Problem.jl")

using Solid
using Fem
using Grid
using Problem
using Output
using Algorithm
using FemMPM
using Material
using Fix
using Basis
using BodyForce

#function main()

    # problem parameters
	density       = 1010e-12
	K             = 121.7         # bulk modulus
	mu            = 26.1          # shear modulus

	E             = 9*K*mu/(3*K+mu)
	nu            = E/2/mu-1.0

	w             = 120.
	l             = 200.

	vel           = 30e3 # mm/s


    # create the grid of a 1 x 1 square, with 20 x 20 cells
    grid  =  Grid2D(0, l, 0, w, 81, 41)
    basis = CPDIQ4Basis()
    #basis = LinearBasis()


    material  = NeoHookeanMaterial(E,nu,density,10000)
    #material  = ElasticMaterial(E,nu,density,0,0)

	c_dil     = sqrt((material.lambda + 2*material.mu)/material.density)
	dt        = grid.dx/c_dil
	dtime     = 0.2 * dt


    solid1   = Solid2D("ring.msh",material)
    solid2   = Solid2D("ring.msh",material)

    v0 = SVector{2,Float64}([vel  0.0])

    # assign initial velocity for the particles
    Solid.assign_velocity(solid1, v0)
    Solid.assign_velocity(solid2,-v0)

    Solid.move_cpdi(solid1,SVector{2,Float64}([ l/4,  w/2]))
    Solid.move_cpdi(solid2,SVector{2,Float64}([ 0.75*l,  0.5*w]))

    solids = [solid1, solid2]

    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("Sound vel : %+.6e \n", c_dil)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))

    Tf      = 0.0047#5e-3 #3.5e-0
    interval= 20

    bodyforce = ConstantBodyForce2D(@SVector[0.,0.])


    data               = Dict()
    data["total_time"] = Tf
    data["dt"]         = dtime
    data["time"]       = 0.
    data["bodyforce"]  = bodyforce


	output2  = OvitoOutput(interval,"rings-cpdi-results/",["vx","sigmaxx"])
	fix      = EnergiesFix(solids,"rings-cpdi-results/energies.txt")
    algo1    = USL(0.)
    algo2    = MUSL(.99)
    
	plotGrid(output2,grid,0)
    plotParticles_2D(output2,solids,[grid.lx, grid.ly],
                         [grid.nodeCountX, grid.nodeCountY],0)
    solve_explicit_dynamics_2D(grid,solids,basis,algo2,output2,fix,data)


    pyFig_RealTime = PyPlot.figure("MPM 2Disk FinalPlot", figsize=(8/2.54, 4/2.54))
	PyPlot.clf()
	pyPlot01 = PyPlot.gca()
	PyPlot.subplots_adjust(left=0.15, bottom=0.25, right=0.65)
	pyPlot01[:grid](b=true, which="both", color="gray", linestyle="-", linewidth=0.5)
	pyPlot01[:set_axisbelow](true)
	#pyPlot01[:set_xlim](0.0, 100.0)
	#pyPlot01[:set_ylim](0.0, 3.0)
	pyPlot01[:set_xlabel]("time (s)", fontsize=8)
	pyPlot01[:set_ylabel]("energy (\$\\times 10^{-3}\$ Nm)", fontsize=8)
	pyPlot01[:set_xticks](collect(0.0:10.0:100.0))
	pyPlot01[:tick_params](axis="both", which="major", labelsize=8)
	#pyPlot01[:set_yticks](collect(0.0:1.0:3.0))
	PyPlot.plot(fix.recordTime, c="blue", fix.kinEnergy, "-", label="\$ K \$", linewidth=1.0)
	#PyPlot.hold(true)
	PyPlot.plot(fix.recordTime, c="red", fix.strEnergy, "-", label="\$ U \$", linewidth=1.0)
	PyPlot.plot(fix.recordTime, c="green", fix.kinEnergy + fix.strEnergy, "-", label="\$ K+U \$", linewidth=1.0)
	PyPlot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=8)
	PyPlot.savefig("plot_2Disk_Julia.pdf")

# end

# @time main()
