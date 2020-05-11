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

function main()
    # problem  material parameters: N, mm and s
	fGravity      = 0.0
	density       = 2400e-12
	youngModulus  = 40000.0
	poissonRatio  = 0.2
	Gf            = 0.11
	l0            = 10.0
    # problem geometry parameters
	L  = 500. # mm
	W  = 100.
	R  = 10.
	r  = 5.
	d  = 5.
	D  = 10.

	vel = 50.  # velocity of the impactor

    numx = 121
    numy = 50
    # create the grid of a 1 x 1 square, with 20 x 20 cells
    grid  = Grid2D(-0.5*L,0.5*L,-0.5*W-2*r,0.5*W+D+R,numx+1, numy+1)
    basis = LinearBasis()

	ppc     = 2
    fOffset = grid.dx/ppc
    dx      = fOffset
    leftCircle   = buildParticleForCircle([-0.5*L+d; -0.5*W-r], r, fOffset)
    rightCircle  = buildParticleForCircle([ 0.5*L-d; -0.5*W-r], r, fOffset)
    impactRec    = buildParticleForRectangle([ 0.; 0.5*W+0.5*R], 4*R, R, fOffset)
    rectangle    = buildParticleForRectangleWithANotch([0.;0.],L, W, fOffset,-0.5*grid.dx,0.5*grid.dx,0.)

    material1    = ElasticMaterial(youngModulus,poissonRatio,density,Gf,l0)
    material2    = RigidMaterial([0.,0.],  density)
    material3    = RigidMaterial([0.,-vel],density)

    solid1       = Solid2D(rectangle,material1)
    solid2       = Solid2D([leftCircle;rightCircle],material2)
    solid3       = Solid2D(impactRec,material3)

    solid1.mass          .*= dx * dx
    solid1.volume        .= dx * dx
    solid1.volumeInitial .= dx * dx

    solid2.mass          .*= dx * dx
    solid2.volume        .= dx * dx
    solid2.volumeInitial .= dx * dx

	solid3.mass          .*= dx * dx
	solid3.volume        .= dx * dx
	solid3.volumeInitial .= dx * dx

    solids = [solid1, solid2, solid3]

	c_dil     = sqrt((material1.lambda + 2*material1.mu)/material1.density)
	dt        = grid.dx/c_dil
	dtime     = 0.5 * dt

    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("grid size : %+.6e \n", grid.dx)
    @printf("time step : %+.6e \n", dtime  )


    Tf       = 3.5 #3.5e-0
    interval = 100

    bodyforce = ConstantBodyForce2D(fGravity)
	output2   = OvitoOutput(interval,"three-point-bending-results/",
	                         ["sigmaxx","damage","crackforce"])
	#fix     =
    problem   = ExplicitDynamicsPhaseField2D(grid,solids,basis,Tf,bodyforce,output2,[])
    algo1     = USL(1e-9)
    algo2     = MUSL(0.99)

	plotParticles(problem.output,solids,[grid.lx, grid.ly],[grid.nodeCountX, grid.nodeCountY],0)
    plotParticles(problem.output,grid,0)
    solve(problem,algo2,dtime)
end

@time main()
