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
    steelRho  = 7850e-12
    steelE    = 200.0e3
    steelNu   = 0.3
    v         = 1160.0e3

    alumRho   = 2700e-12
    alumE     = 78.2e3
    alumNu    = 0.3
    alumFy    = 300.
    alumK     = 0.    # hardening modulus


	# grid creation
	grid  = Grid2D(0.,60.0, 0.,60.0, 61, 61)
	basis = LinearBasis()

    solid1   = FEM2D("disk-steel.msh")
    solid2   = FEM2D("rect.msh")


    material1 = ElasticMaterial(steelE,steelNu,steelRho,0.,0.)
    material2 = ElastoPlasticMaterial(alumE,alumNu,alumRho,alumFy,alumK,solid2.parCount)

    v0 = SVector{2,Float64}([0.0 -v])

    # assign initial velocity for the particles
    Fem.assign_velocity(solid1, v0)
    Fem.move(solid1,SVector{2,Float64}([ 30.,  40. + 10.]))
   # Fem.move(solid2,SVector{2,Float64}([ 1.,  1.]))
    

    solids = [solid1 solid2]
    mats   = [material1 material2]

    fixNodes(solid2,"bottom")
    fixNodes(solid2,"left-right")

    fixXForLeft(grid)
    fixYForLeft(grid)
    fixXForRight(grid)
    fixYForRight(grid)
    fixXForBottom(grid)
    fixYForBottom(grid)

    Tf      = 50e-6
    interval= 200
	dtime   = 1.0e-8

    # c_dil     = sqrt(alumE/alumRho)
    # dt        = grid.dy/c_dil
    # dtime     = 0.1 * dt


	output2  = VTKOutput(interval,"impact-femp-results/",["vx","sigmaxx"])
	fix      = DisplacementFemFix(solid2,"impact-femp-results/",212)
    algo1    = USL(0.)
    algo2    = TLFEM(0.,0.999)
    bodyforce = ConstantBodyForce2D([0.,0.])



    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))

	plotGrid(output2,grid)
	plotParticles_2D(output2,solids,mats,0)
    solve_explicit_dynamics_femp_2D(grid,solids,mats,basis,bodyforce,algo2,output2,fix,Tf,dtime)


# end

# @time main()
