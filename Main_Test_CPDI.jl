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


push!(LOAD_PATH,"./")

# import Gadfly
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")
using DelimitedFiles



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

#function main()
	fGravity  = 1000e3

	steelRho  = 1050e-12
    steelE    = 1.
    steelNu   = 0.3

	c         = sqrt(steelE/steelRho)

	# grid creation and basis
	grid  = Grid2D(0.,2000., 0., 3500., 5, 8)
	basis = CPDIQ4Basis()

    # do not forget to shift the solid to the right one cell (to have ghost cell)
	nodes,elems = createMeshForRectangle([0. 0.],[1000. 0.],
	                                     [1000. 1000.],[0. 1000.],6,6)
    nodes[1,:] .+= 500 #grid.dx
    nodes[2,:] .+= 3500-1000 - grid.dy

	material    = NeoHookeanMaterial(steelE,steelNu,steelRho,length(nodes))
	#println(nodes)
	solid1      = Solid2D(nodes,elems,material)

    solids = [solid1]

	@printf("	Disk,   number of material points: %d \n", solid1.parCount)

	bodyforce = ConstantBodyForce2D(@SVector[0.,-fGravity])
 
    algo1    = USL(1e-9)
    algo2    = MUSL(1.)


    dtime   = 0.1*grid.dx/c;
    Tf      = 0.25
    interval= 1


    data               = Dict()
    data["total_time"] = Tf
    data["dt"]         = dtime
    data["time"]       = 0.
    data["dirichlet_grid"] = [("top",(0,1))] # => fix left edge of the grid
    data["bodyforce"]  = bodyforce
    data["ghostcell"]  = true

    dir      = "cpdi-test-results/"

	output2  = OvitoOutput(interval,dir,["sigmaxx"])
	fix      = DisplacementFix(solids,@SVector[1083.3333333333333, 2250.0],dir)

 
    algo1    = USL(1e-9)
    algo2    = MUSL(1.)

	report(grid,solids,dtime)
	plotGrid(output2,grid,0)

	#plotParticles(problem.output,solids,[grid.lx, grid.ly],[grid.nodeCountX, grid.nodeCountY],0)
    #plotParticles_2D(output2,grid,0)

	#reset_timer!()
    solve_explicit_dynamics_2D(grid,solids,basis,algo2,output2,fix,data)

    v = readdlm("cpdi-test-results/recorded-position.txt")

    PyPlot.plot(v[:, 1], v[:, 2])
    PyPlot.plot(v[:, 1], v[:, 3])

# end

# @time main()
