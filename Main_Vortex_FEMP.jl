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

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/juMP")
#
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")


#include("./Grid.jl")
#include("./Problem.jl")

using Solid
using Fem
using Grid
using Problem
using Output
using Algorithm
using FemMPM
using Material
using Fix
using Basis
using BodyForce

#function main()

    # problem parameters
	density       = 1050.
	E             = 1e6
	nu            = 0.3
	G             = 0.2
	T             = 1.

	mu     = E/(2*(1+nu))
    lambda = E*nu/((1+nu)*(1-2*nu))

    Ro = 1.25
    Ri = 0.75
    l  = 2.8/2


    # create the grid of a 1 x 1 square, with 20 x 20 cells
    #grid  =  Grid2D(-l, l, -l, l, 81, 81)
    grid  =  Grid2D(0, 2*l, 0, 2*l, 101, 81)
    basis = QuadBsplineBasis()
    #basis = LinearBasis()


    material  = NeoHookeanMaterial(E,nu,density)
    #material  = ElasticMaterial(E,nu,density,0,0)

	c_dil     = sqrt((material.lambda + 2*material.mu)/material.density)
	dt        = grid.dx/c_dil
	dtime     = 0.05 * dt


    #solid1   = FEM2D("vortex.msh",material)
    solid1   = FEM2D("vortex-regular.msh",material)
    

    # v0 = SVector{2,Float64}([vel  0.0])

    # # assign initial velocity for the particles
    # Fem.assign_velocity(solid1, v0)
    # Fem.assign_velocity(solid2,-v0)

    Fem.move(solid1,SVector{2,Float64}([ l,  l]))
    
    fixNodes(solid1, "inner")
    fixNodes(solid1, "outer")

    solids = [solid1]

    bodyforce = VortexBodyForce2D(G,T,Ri,Ro,density,mu)

    @printf("Total number of material points: %d \n", solid1.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("Sound vel : %+.6e \n", c_dil)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))

    Tf      = 1.#5e-3 #3.5e-0
    interval= 100

	output2  = VTKOutput(interval,"vortex-femp-results/",["vx","sigmaxx"])
	fix      = EnergiesFix(solids,"vortex-femp-results/energies.txt")
    
    algo1    = TLFEM(0.)
    algo2    = TLFEMFull(1e-8)
    
	plotGrid(output2,grid)
	plotParticles_2D(output2,solids,0)
    solve_explicit_dynamics_femp_2D(grid,solids,basis,bodyforce,algo1,output2,fix,Tf,dtime)


 

# end

# @time main()
