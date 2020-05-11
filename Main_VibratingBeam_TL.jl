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
using Basis
using Problem
using Output
using Fix
using BodyForce

function main()
	fGravity  = 10.0

	steelRho  = 1050.
    steelE    = 1e6
    steelNu   = 0.3

	c         = sqrt(steelE/steelRho)

	# grid creation
	grid1 = Grid2D(0.,4.0, 0.,1.0, 41, 11)
	basis = LinearBasis()

    # how to calculate offset: cell size = 60/50, ppc = 3
	fOffset = grid1.dx/3.0
	coords  = buildParticleForRectangle([2.0; .5], 4.0, 1., fOffset)

	material = NeoHookeanMaterial(steelE,steelNu,steelRho)

	solid1   = Solid2D(coords,material)


    solid1.mass          .*=fOffset * fOffset
    solid1.volume        .= fOffset * fOffset
    solid1.volumeInitial .= fOffset * fOffset

    solids = [solid1]
	grids  = [grid1]

	@printf("	Disk,   number of material points: %d \n", solid1.parCount)

    # Boundary conditions
    fixXForLeft(grid1)
    fixYForLeft(grid1)

    dtime   = 0.2*grid1.dx/c;
    Tf      = 3.
    interval= 100

	bodyforce = ConstantBodyForce2D(fGravity)

	#output1  = PyPlotOutput(interval,"impact-results/","Impact",(4., 4.))

	# output2  = OvitoOutput(interval,"vibrating-beam-TL-FEM-results/",["sigmaxx", "sigmayy"])
	# fix      = DisplacementFix(solids,[3.98333 0.0166625],"vibrating-beam-TL-FEM-results/")
    # algo     =  FEM(1e-9)
	output2  = OvitoOutput(interval,"vibrating-beam-TL-results/",["sigmaxx", "sigmayy"])
	fix      = DisplacementFix(solids,[3.98333 0.0166625],"vibrating-beam-TL-results/")
    algo     =  USL(1e-9)
    problem  = ExplicitDynamicsTL2D(grids,solids,basis,Tf,bodyforce,output2,[fix])


    solve(problem, algo, dtime)
end

@time main()
