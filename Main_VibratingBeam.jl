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

function main()
	fGravity  = 10.0

	steelRho  = 1050.
    steelE    = 10.e6
    steelNu   = 0.3

	c         = sqrt(steelE/steelRho)

	# grid creation
	grid = Grid2D(5.0, 10.0, 501, 101)

    # how to calculate offset: cell size = 60/50, ppc = 3
	fOffset = 5.0/100/3.0
	coords  = buildParticleForRectangle([2.0; 5.0], 4.0, 1., fOffset)

	material = NeoHookeanMaterial(steelE,steelNu,steelRho)

	solid1   = Solid2D(coords,material)


    solid1.mass          .*=fOffset * fOffset
    solid1.volume        .= fOffset * fOffset
    solid1.volumeInitial .= fOffset * fOffset

    solids = [solid1]

	@printf("	Disk,   number of material points: %d \n", solid1.parCount)

    # Boundary conditions
    fixXForLeft(grid)
    fixYForLeft(grid)

    dtime   = 0.2*grid.dx/c;
    Tf      = 3.
    interval= 200

	bodyforce = ConstantBodyForce2D(fGravity)

	#output1  = PyPlotOutput(interval,"impact-results/","Impact",(4., 4.))
	output2  = OvitoOutput(interval,"vibrating-beam-results/",[])
    problem  = ExplicitSolidMechanics2D(grid,solids,Tf,bodyforce,output2,[])
    algo     =  MUSL()
    solve(problem, algo, dtime)
end

@time main()
