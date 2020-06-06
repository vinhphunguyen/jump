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
	density       = 8940e-12
	E             = 115e3
	nu            = 0.31
	A             = 65.
	B             = 356.
	C             = 0.013
	n             = 0.37
	eps0dot       = 1.0

	vel           = 190e3 # mm/s


    # create the grid of a 1 x 1 square, with 20 x 20 cells
    grid  =  Grid3D(0, 14, 0, 28, 0, 14, 20, 60, 20)
    basis = QuadBsplineBasis()
    #basis = LinearBasis()

    # (E,nu,density,A,B,C,n,eps0dot,cellsize,parCount)
    #material  = JohnsonCookMaterial(E,nu,density,A,B,C,n,eps0dot,.42,14160)
    #material  = NeoHookeanMaterial(E,nu,density)
    material  = ElastoPlasticMaterial(E,nu,density,A,B,14160)
    
	c_dil     = sqrt(E/density)
	dt        = grid.dy/c_dil
	dtime     = 0.1 * dt


    solid1   = FEM3D("taylor-bar.msh",material)

    v0       = SVector{3,Float64}([0.0 -vel 0.])

    # assign initial velocity for the particles
    Fem.assign_velocity(solid1, v0)
    
    Fem.move(solid1,SVector{3,Float64}([ 7, 1.6, 7]))
    

    solids = [solid1]

   
 
    fixYForBottom(grid)
    

    @printf("Total number of material points: %d \n", solid1.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("Sound vel : %+.6e \n", c_dil)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))

    Tf      = 63e-3
    interval= 500

	output   = VTKOutput(interval,"taylor-femp-results/",["vx","sigmaxx"])
	fix      = EnergiesFix(solids,"taylor-femp-results/energies.txt")
    algo1    = USL(0.)
    algo2    = TLFEM(0.)
    bodyforce = ConstantBodyForce3D([0., 0.,0.])

	plotGrid(output,grid,0)
	plotParticles_3D(output,solids,0)
    solve_explicit_dynamics_femp_3D(grid,solids,basis,bodyforce,algo2,output,fix,Tf,dtime)

    D = abs(solid1.pos[2][1] - solid1.pos[4][1])
    L = abs(solid1.pos[2][2] - solid1.pos[6][2])
    #W = abs(solid1.pos[?][1] - solid1.pos[?][1]) # 5.08 from the bottom 
# end

# @time main()
