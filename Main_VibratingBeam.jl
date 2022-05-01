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
# Solved with the ULMPM
# Output in folder "vibratingbeam-mpm-results/", with lammps dump files and energies.txt

push!(LOAD_PATH,"./")


# import Gadfly
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")


#include("./Grid.jl")
#include("./Problem.jl")

using Material
using Algorithm
using Solid
using Grid
using Problem
using Output
using Fix
using BodyForce
using Basis
using Util

function main()
	fGravity  = 10.0

	steelRho  = 1050.
    steelE    = 1e6
    steelNu   = 0.3

	c         = sqrt(steelE/steelRho)

	# grid creation
	grid   = Grid2D(0.,5.0, 0.,10.0, 51, 101)
	basis1 = LinearBasis()
    basis2 = QuadBsplineBasis()
    basis3 = CubicBsplineBasis()

    # how to calculate offset: cell size = 60/50, ppc = 3
	fOffset = 5.0/50/3.0
	coords  = buildParticleForRectangle([2.0; 5.0], 4.0, 1., fOffset)

	material = NeoHookeanMaterial(steelE,steelNu,steelRho,length(coords))

	solid1   = Solid2D(coords,material)

    # solid.mass already computed from density of material, now need volume
    solid1.mass          .*=fOffset * fOffset 
    solid1.volume        .= fOffset * fOffset
    solid1.volumeInitial .= fOffset * fOffset

    solids = [solid1]

	@printf("	Disk,   number of material points: %d \n", solid1.parCount)

    dtime   = 0.4*grid.dx/c;
    Tf      = 3.
    interval= 1

	bodyforce = ConstantBodyForce2D(@SVector[0.,-fGravity])

    data               = Dict()
    data["total_time"] = Tf
    data["dt"]         = dtime
    data["time"]       = 0.
    data["dirichlet_grid"] = [("left",(1,1))] # => fix left edge of the grid
    data["bodyforce"]  = bodyforce


	#output1  = PyPlotOutput(interval,"impact-results/","Impact",(4., 4.))
	output2  = OvitoOutput(interval,"vibratingbeam-mpm-results/",["sigmaxx"])
	fix      = EmptyFix()#DisplacementFix(solids,@SVector[3.9916666666666667,4.508333333333334],2)

 
    algo1    = USL(1e-9)
    algo2    = MUSL(1.)

	report(grid,solids,dtime)
	plotGrid(output2,grid,0)

	#plotParticles(problem.output,solids,[grid.lx, grid.ly],[grid.nodeCountX, grid.nodeCountY],0)
    #plotParticles_2D(output2,grid,0)

	#reset_timer!()
    solve_explicit_dynamics_2D(grid,solids,basis3,algo1,output2,fix,data)
end

@time main()
