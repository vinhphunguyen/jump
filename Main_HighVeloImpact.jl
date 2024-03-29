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

# Input file for the vibrating cantilever beam proposed by Brannon et al.
# Solved with the CPDI-Q4
# Output in folder "vibratingbeam-cpdi-results/", with lammps dump files and energies.txt

push!(LOAD_PATH,"./")


# import Gadfly
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")
using BenchmarkTools
using TimerOutputs

#include("./Grid.jl")
#include("./Problem.jl")


using Solid
using Grid
using Problem
using Output
using Fix
using Basis
using BodyForce
using Material
using Algorithm
using Util

function main()
	fGravity  = 0.0

	steelRho  = 7850e-12
    steelE    = 200.0e3
    steelNu   = 0.3
    v         = 1160.0e3

    alumRho   = 2700e-12
    alumE     = 78.2e3
    alumNu    = 0.3
    alumFy    = 300.
    alumK     = 0.    # hardening modulus


	# grid creation
	grid  = Grid2D(0.,60.0, 0.,60.0, 51, 51)
	basis = LinearBasis()

    # how to calculate offset: cell size = 60/50, ppc = 3
	fOffset = 60.0/50/3.0
	coords1 = buildParticleForCircle([30.0; 50.0], 9.6/2.0, fOffset)
	coords2 = buildParticleForRectangle([30.0; 20.0], 60.0, 40., fOffset)

	material1 = ElasticMaterial(steelE,steelNu,steelRho,0.,0.)
	material2 = ElastoPlasticMaterial(alumE,alumNu,alumRho,alumFy,alumK,length(coords2))


	solid1 = Solid2D(coords1,material1)
    solid2 = Solid2D(coords2,material2)

    solid1.mass          .*=fOffset * fOffset
    solid1.volume        .= fOffset * fOffset
    solid1.volumeInitial .= fOffset * fOffset

    solid2.mass          .*=fOffset * fOffset
    solid2.volume        .= fOffset * fOffset
    solid2.volumeInitial .= fOffset * fOffset

    v0 = SVector{2,Float64}([ 0 -v])

    # assign initial velocity for the particles
    assign_velocity(solid1, v0)


    solids = [solid1, solid2]

	bodyforce = ConstantBodyForce2D(fGravity)


    # Boundary conditions

    fixXForBottom(grid)
    fixYForBottom(grid)
    fixXForLeft(grid)
    fixYForLeft(grid)
    fixXForRight(grid)
    fixYForRight(grid)


    Tf      = 48.0e-6
    interval= 200
	dtime   = 1.0e-8

	output1  = PyPlotOutput(interval,"impact-results/","Impact",(4., 4.))
	output2  = OvitoOutput(interval,"impact-results/",["pstrain", "vonMises"])
	fix      = DisplacementFix(solids,[29.8 40.],"impact-results/")
    algo1     =  MUSL(1.)
    algo2     =  USL(1e-11)

	report(grid,solids,dtime)

	# reset_timer!()
	# @time solve_explicit_dynamics_2D(grid,solids,basis,algo2,output2,fix,Tf,dtime)
    # print_timer()

	solve_explicit_dynamics_2D(grid,solids,basis,algo1,output2,fix,Tf,dtime)
end

@time main()
