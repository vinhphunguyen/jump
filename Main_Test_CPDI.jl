# Phu Nguyen, Monash University
# 20 March, 2020 (Coronavirus outbreak)

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/jMPM/src")
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
using Mesh

#function main()
	fGravity  = 1000e3

	steelRho  = 1050e-12
    steelE    = 1.
    steelNu   = 0.3

	c         = sqrt(steelE/steelRho)

	# grid creation and basis
	grid  = Grid2D(2000., 3500., 5, 8)
	basis = CPDIQ4Basis()

    # do not forget to shift the solid to the right one cell (to have ghost cell)
	nodes,elems = createMeshForRectangle([0. 0.],[1000. 0.],
	                                     [1000. 1000.],[0. 1000.],6,6)
    nodes[1,:] .+= grid.dx
    nodes[2,:] .+= 4*grid.dx

	material    = NeoHookeanMaterial(steelE,steelNu,steelRho)
	#println(nodes)
	solid1      = Solid2D(nodes,elems,material)

    solids = [solid1]

	@printf("	Disk,   number of material points: %d \n", solid1.parCount)

    # Boundary conditions
    fixXForTop(grid,ghostcell=true)
    fixYForTop(grid,ghostcell=true)

    dtime   = 0.1*grid.dx/c;
    Tf      = 0.25
    interval= 10

	bodyforce = ConstantBodyForce2D(fGravity)

	#output1  = PyPlotOutput(interval,"impact-results/","Impact",(4., 4.))
	output2  = OvitoOutput(interval,"cpdi-test-results/",[])
    problem  = ExplicitSolidMechanics2D(grid,solids,basis,Tf,bodyforce,output2,[])
    algo     =  MUSL()
    solve(problem, algo, dtime)
# end
#
# @time main()
