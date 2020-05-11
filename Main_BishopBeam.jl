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

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/jMPM/src")
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

#Simulating the pervasive fracture of materials and structures using randomly close packed Voronoi tessellations
#Joseph E. Bishop

function main()
    # problem  material parameters: N, mm and s
	fGravity      = 0.0
	density       = 2250e-12
	youngModulus  = 28.3e3
	poissonRatio  = 0.2
	Gf            = 0.005#0.057
	l0            = 20.0
    # problem geometry parameters
	L  = 1.83e3 # mm
	W  = 0.3e3


	vel  = 7.6e3  # velocity of the impactor

    numx = 200
    numy = 50
    # create the grid of a 1 x 1 square, with 20 x 20 cells
    grid  = Grid2D(0,2*L,0,1.1*L,numx+1, numy+1)
    basis = LinearBasis()

	ppc     = 2
    fOffset = grid.dx/ppc
    dx      = fOffset

    rectangle1   = buildParticleForRectangle([0.5*L;0.5*W],L, W, fOffset)
    rectangle2   = buildParticleForRectangle([1.0*L;0.5*dx],1.5*L, dx, fOffset)

    material1    = ElasticMaterial(youngModulus,poissonRatio,density,Gf,l0)
    material2    = RigidMaterial([0.,0.],  density)


    solid1       = Solid2D(rectangle1,material1)
    solid2       = Solid2D(rectangle2,material2)

    rotate(solid1,-45.)
    move(solid1,[1000.,0.6*L])

	assign_velocity(solid1,[0.,-vel])

    solid1.mass          .*= dx * dx
    solid1.volume        .= dx * dx
    solid1.volumeInitial .= dx * dx

    solid2.mass          .*= dx * dx
    solid2.volume        .= dx * dx
    solid2.volumeInitial .= dx * dx


    solids = [solid1, solid2]

    Tf       = 3.5 #3.5e-0
    interval = 100

	c_dil     = sqrt((material1.lambda + 2*material1.mu)/material1.density)
	dt        = grid.dx/c_dil
	dtime     = 0.2 * dt


    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("dtime : %+.6e \n", dtime)
    @printf("dx    : %+.6e \n", grid.dx)



    bodyforce = ConstantBodyForce2D(fGravity)
	output2   = OvitoOutput(interval,"Bishop-beam-results/",
	                         ["sigmaxx","damage","crackforce"])
	#fix     =
    problem   = ExplicitDynamicsPhaseField2D(grid,solids,basis,Tf,bodyforce,output2,[])
    algo1     = USL(1e-9)
    algo2     = MUSL(0.999999)

	plotParticles(problem.output,solids,[grid.lx, grid.ly],[grid.nodeCountX, grid.nodeCountY],0)
    plotParticles(problem.output,grid,0)
    solve(problem,algo2,dtime)
end

@time main()
