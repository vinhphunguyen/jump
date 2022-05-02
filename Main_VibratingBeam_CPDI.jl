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
using DelimitedFiles




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
using Mesh
using Util

function main()
	fGravity  = 10.0

	steelRho  = 1050.
    steelE    = 1e6
    steelNu   = 0.3

	c         = sqrt(steelE/steelRho)

	# grid creation and basis
	grid  = Grid2D(0,5.0, 0, 10.0, 17, 34)
	basis = CPDIQ4Basis()

    # do not forget to shift the solid to the right one cell (to have ghost cell)
	nodes,elems = createMeshForRectangle([1.1*grid.dx 4.5],[4.0+1.1*grid.dx 4.5],
	                                     [4.0+1.1*grid.dx 5.5],[1.1*grid.dx 5.5],40,10)
	material    = NeoHookeanMaterial(steelE,steelNu,steelRho,size(elems,1))
	#println(nodes)
	solid1      = Solid2D(nodes,elems,material)

    solids = [solid1]

	@printf("	Disk,   number of material points: %d \n", solid1.parCount)

    dtime   = 0.2*grid.dx/c;
    Tf      = 3.
    interval= 1

	bodyforce = ConstantBodyForce2D(@SVector[0.,-fGravity])

    data               = Dict()
    data["total_time"] = Tf
    data["dt"]         = dtime
    data["time"]       = 0.
    data["dirichlet_grid"] = [("left",(1,1))] # => fix left edge of the grid
    data["bodyforce"]  = bodyforce
    data["ghostcell"]  = true

	output2  = OvitoOutput(interval,"vibratingbeam-cpdi-results/",["sigmaxx"])
	fix      = EmptyFix()#DisplacementFix(solids,@SVector[4.060000000000000,4.55],2)

 
    algo1    = USL(1e-12)
    algo2    = MUSL(1.)

	report(grid,solids,dtime)
	plotGrid(output2,grid,0)

	#plotParticles(problem.output,solids,[grid.lx, grid.ly],[grid.nodeCountX, grid.nodeCountY],0)
    #plotParticles_2D(output2,grid,0)

	#reset_timer!()
    solve_explicit_dynamics_2D(grid,solids,basis,algo1,output2,fix,data)

    v = readdlm("vibratingbeam-cpdi-results/recorded-position.txt")

    PyPlot.plot(v[:, 1], v[:, 2])
    PyPlot.plot(v[:, 1], v[:, 3])
end
#
@time main()
