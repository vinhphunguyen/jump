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
    steelRho  = 7850e-9
    steelE    = 200.0
    steelNu   = 0.3
    v         = 1160.0

    alumRho   = 2700e-9
    alumE     = 78.2
    alumNu    = 0.3
    alumFy    = 300.
    alumK     = 0.    # hardening modulus


    A       = .3*1.224744871391589
    B       = 0.0
    C       = 0#0.013
    n       = 0.37
    eps0dot = 1e-3
    Tm      = 1600.


	# grid creation
	grid  = Grid3D(0.,61.0, 0.,60.0, 0,61, 71, 61, 71)
	basis = LinearBasis()

    solid1   = FEM3D("sphere.msh")
    solid2   = FEM3D("rect-3D.msh")
    #solid2   = FEM3D("cylinder.msh")
    #solid2   = FEM3D("cylinder.inp")


    material1 = ElasticMaterial(steelE,steelNu,steelRho,0.,0.)
    material2 = JohnsonCookMaterial(alumE,alumNu,alumRho,A,B,C,n,eps0dot,.42,solid2.parCount) 
    #ElastoPlasticMaterial(alumE,alumNu,alumRho,alumFy,alumK,solid2.parCount)

    v0 = SVector{3,Float64}([0.0 -v 0.])

    # assign initial velocity for the particles
    Fem.assign_velocity(solid1, v0)
    Fem.move(solid1,SVector{3,Float64}([ 30.,  40. + 10., 30.]))
    #Fem.move(solid2,SVector{3,Float64}([ 35.,  0., 35.]))
    Fem.move(solid2,SVector{3,Float64}([ .5,  0., .5]))
    

    solids = [solid1 solid2]
    mats   = [material1 material2]

    #fixNodes(solid2,"bottom")
    #fixNodes(solid2,"left-right")

    fixXForLeft(grid)
    fixYForLeft(grid)
    fixZForLeft(grid)
    fixXForRight(grid)
    fixYForRight(grid)
    fixZForRight(grid)
    fixXForBottom(grid)
    fixYForBottom(grid)
    fixZForBottom(grid)
    fixXForFront(grid)
    fixYForFront(grid)
    fixZForFront(grid)
    fixXForBack(grid)
    fixYForBack(grid)
    fixZForBack(grid)

    #fixForTop(grid)
  


    c_dil     = sqrt(alumE/alumRho)
    dt        = grid.dy/c_dil
    dtime     = 0.1 * dt


    Tf      = 50e-3
    interval= 200
    #dtime   = 1.0e-8


	output2  = VTKOutput(interval,"impact-femp-3D-results/",["vx","sigmaxx"])
	fix      = DisplacementFemFix(solid2,"impact-femp-3D-results/",2)
    algo1    = USL(0.)
    algo2    = TLFEM(0.,0.999)
    bodyforce = ConstantBodyForce3D([0.,0.,0.])



    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))

	plotGrid(output2,grid,0)
	plotParticles_3D(output2,solids,mats,0)
    solve_explicit_dynamics_femp_3D(grid,solids,mats,basis,bodyforce,algo2,output2,fix,Tf,dtime)


    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))
# end

# @time main()
