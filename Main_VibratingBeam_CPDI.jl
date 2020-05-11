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
	fGravity  = 10.0

	steelRho  = 1050.
    steelE    = 10.e6
    steelNu   = 0.3

	c         = sqrt(steelE/steelRho)

	# grid creation and basis
	grid  = Grid2D(5.0, 10.0, 51, 101)
	basis = CPDIQ4Basis()

    # do not forget to shift the solid to the right one cell (to have ghost cell)
	nodes,elems = createMeshForRectangle([1.1*grid.dx 4.5],[4.0+1.1*grid.dx 4.5],
	                                     [4.0+1.1*grid.dx 5.5],[1.1*grid.dx 5.5],40,10)
	material    = NeoHookeanMaterial(steelE,steelNu,steelRho)
	#println(nodes)
	solid1      = Solid2D(nodes,elems,material)

    solids = [solid1]

	@printf("	Disk,   number of material points: %d \n", solid1.parCount)

    # Boundary conditions
    fixXForLeft(grid,ghostcell=true)
    fixYForLeft(grid,ghostcell=true)

    dtime   = 0.2*grid.dx/c;
    Tf      = 3.
    interval= 200

	bodyforce = ConstantBodyForce2D(fGravity)

	#output1  = PyPlotOutput(interval,"impact-results/","Impact",(4., 4.))
	output2  = OvitoOutput(interval,"vibratingbeam-cpdi-results/",[])
    problem  = ExplicitSolidMechanics2D(grid,solids,basis,Tf,bodyforce,output2,[])
    algo     =  MUSL()
    solve(problem, algo, dtime)
# end
#
# @time main()
