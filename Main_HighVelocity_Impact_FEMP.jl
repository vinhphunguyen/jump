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
	grid  = Grid2D(0.,62.0, 0.,62.0, 31, 31)
	basis = LinearBasis()


	material1 = ElasticMaterial(steelE,steelNu,steelRho,0.,0.)
	material2 = ElastoPlasticMaterial(alumE,alumNu,alumRho,alumFy,alumK,9600)

    solid1   = FEM2D("disk-steel.msh",material1)
    solid2   = FEM2D("rect.msh",material2)

    v0 = SVector{2,Float64}([0.0 -v])

    # assign initial velocity for the particles
    Fem.assign_velocity(solid1, v0)
    Fem.move(solid1,SVector{2,Float64}([ 30.,  40. + 10.]))
    Fem.move(solid2,SVector{2,Float64}([ 1.,  1.]))
    

    solids = [solid1 solid2]

    fixNodes(solid2,"bottom")
    fixNodes(solid2,"left-right")

     fixXForLeft(grid)
    fixYForLeft(grid)
    fixXForRight(grid)
    fixYForRight(grid)


    Tf      = 48.0e-6
    interval= 20
	dtime   = 5.0e-9

	output2  = VTKOutput(interval,"impact-femp-results/",["vx","sigmaxx"])
	fix      = EnergiesFix(solids,"impact-femp-results/energies.txt")
    algo1    = USL(0.)
    algo2    = TLFEM(0.)
    bodyforce = ConstantBodyForce2D([0.,0.])



    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))

	plotGrid(output2,grid)
	plotParticles_2D(output2,solids,0)
    solve_explicit_dynamics_femp_2D(grid,solids,basis,bodyforce,algo2,output2,fix,Tf,dtime)


# end

# @time main()
