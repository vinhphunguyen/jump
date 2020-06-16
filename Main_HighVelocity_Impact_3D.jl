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
    v         = 2000.0

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
	grid  = Grid3D(0.,60.0, 0.,60.0, 0,60, 71, 61, 71)
	basis = LinearBasis()

    dx1 = fOffset1 = grid.dx/2
    dx2 = fOffset2 = grid.dx/1

    coords1 = buildParticleForSphere([0.0; 0.0; 0.0], 9.53/2, fOffset1)
    #coords2 = buildParticleForBlock([0.0; 0.0; 0.0], 60., 40., 60., fOffset2)
    coords2 = buildParticleForCylinder([0.;  0.; 0.], 30., 20., fOffset2, fOffset2)
    
    #solid2   = FEM3D("cylinder.msh")
    #solid2   = FEM3D("cylinder.inp")


    material1 = ElasticMaterial(steelE,steelNu,steelRho,0.,0.)
    material2 = JohnsonCookMaterial(alumE,alumNu,alumRho,A,B,C,n,eps0dot,.42,length(coords2)) 
    #ElastoPlasticMaterial(alumE,alumNu,alumRho,alumFy,alumK,solid2.parCount)

    solid1   = Solid3D(coords1,material1)
    solid2   = Solid3D(coords2,material2)

    solid1.mass          .*= dx1^3
    solid1.volume        .=  dx1^3
    solid1.volumeInitial .=  dx1^3

    solid2.mass          .*= dx2^3
    solid2.volume        .=  dx2^3
    solid2.volumeInitial .=  dx2^3

    v0 = SVector{3,Float64}([0.0 -v 0.])

    # assign initial velocity for the particles
    Solid.assign_velocity(solid1, v0)
    Solid.move(solid1,SVector{3,Float64}([ 30.,  40. + 10., 30.]))
    Solid.move(solid2,SVector{3,Float64}([ 30.,  20., 30.]))
    

    solids = [solid1, solid2]
    

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


	output2  = OvitoOutput(interval,"impact-3D-results/",["pressure","vonMises"])
	fix      = DisplacementFix(solids,@SVector[30.42857142857142, 39, 30.42857142857142],"impact-3D-results/")
    algo1    = USL(0.)
    algo2    = MUSL(0.999)
    bodyforce = ConstantBodyForce3D([0.,0.,0.])



    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))

	plotGrid(output2,grid,0)
	plotParticles_3D(output2,solids,[grid.lx, grid.ly, grid.lz],
                    [grid.nodeCountX, grid.nodeCountY, grid.nodeCount],0)
    
    solve_explicit_dynamics_3D(grid,solids,basis,algo2,output2,fix,Tf,dtime)


    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))
# end

# @time main()
