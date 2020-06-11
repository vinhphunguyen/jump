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
	density       = 1010e-9
	K             = 121.7e-3         # bulk modulus
	mu            = 26.1e-3          # shear modulus

	E             = 9*K*mu/(3*K+mu)
	nu            = E/2/mu-1.0

	w             = 200.
	l             = 200.

	vel           = -30 # mm/ms


    # create the grid of a 1 x 1 square, with 20 x 20 cells
    grid  =  Grid2D(0, l, 0, w, 81, 81)
    #basis = QuadBsplineBasis()
    basis = LinearBasis()


    material1  = NeoHookeanMaterial(E,nu,density)
    material2 = RigidMaterial(0,vel,E)

	c_dil     = sqrt((material1.lambda + 2*material1.mu)/material1.density)
	dt        = grid.dx/c_dil
	dtime     = 0.2 * dt


    solid1   = FEM2D("ring.msh",material1)
    solid2   = FEM2D("ring.msh",material1)
    solid3   = FEM2D("ring.msh",material1)
    solid4   = FEM2D("ring.msh",material1)
    solid5   = FEM2D("platen.msh",material2)

    v0 = SVector{2,Float64}([vel  0.0])

    # assign initial velocity for the particles
    Fem.assign_velocity(solid5, v0)
  

    Fem.move(solid1,SVector{2,Float64}([ 60,  40]))
    Fem.move(solid2,SVector{2,Float64}([ 60+80.,  40]))
    Fem.move(solid3,SVector{2,Float64}([ 60,  40+80]))
    Fem.move(solid4,SVector{2,Float64}([ 60+80,  40+80]))
    Fem.move(solid5,SVector{2,Float64}([ 0,  166]))

    solids = [solid1, solid2, solid3, solid4, solid5]

    fixXForLeft(grid)
    fixXForRight(grid)
    fixYForBottom(grid)

    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("Sound vel : %+.6e \n", c_dil)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))

    Tf      = 3.5#5e-3 #3.5e-0
    interval= 10

	output2  = VTKOutput(interval,"four-rings-femp-results/",["vx","sigmaxx"])
	fix      = EnergiesFix(solids,"four-rings-femp-results/energies.txt")
    algo1    = USL(0.)
    algo2    = TLFEM(0.)
    bodyforce = ConstantBodyForce2D([0.,0.])

	plotGrid(output2,grid)
	plotParticles_2D(output2,solids,0)
    solve_explicit_dynamics_femp_2D(grid,solids,basis,bodyforce,algo2,output2,fix,Tf,dtime)



# end

# @time main()
