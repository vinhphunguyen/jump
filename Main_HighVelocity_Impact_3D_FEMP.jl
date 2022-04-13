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


# This file is for the high velocity impact of a sphere into a cylinder sample presented in
# paper 'A generalized particle in cell metghod for explicit solid dynamics', CMAME, V.P. Nguyen et all 2020.

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
    m       = 0.
    χ       = 0.
    Cp      = 0.


	# grid creation: 60x60x60 with origin at 0,0,0 by 70x60x70 cells
    # the actual grid length in x-dir is 60.01 so that there is no FE nodes exactly locating on the grid line!!!
    # similar for z-dir
    # See Basis.jl, getNearestGridPoints(points,start,xp,grid::Grid2D), floor(Int64,(xp[1]-grid.xmin) * fLength_Cell_x + 1.)
    # linear basis for the grid
	grid  = Grid3D(0.,60.01, 0.,60.0, 0.,60.01, 71, 61, 71)
	basis = LinearBasis()

    # solids from meshes: solid1 = sphere which is Gmsh and 
    # solid2=cylinder created in Abaqus (inp file) or Gmsh (msh file)
    solid1   = FEM3D("sphere-impact.msh")
    #solid2   = FEM3D("cylinder.inp")
    solid2   = FEM3D("cylinder-impact.msh")

    # materials: 2, one for each solid
    material1 = ElasticMaterial(steelE,steelNu,steelRho,0.,0.,solid1.parCount)
    material2 = JohnsonCookMaterial(alumE,alumNu,alumRho,A,B,C,n,eps0dot,m,χ,Cp,.42,solid2.parCount) 
    #ElastoPlasticMaterial(alumE,alumNu,alumRho,alumFy,alumK,solid2.parCount)

    v0 = SVector{3,Float64}([0.0 -v 0.])

    # assign initial velocity for the particles
    Fem.assign_velocity(solid1, v0)

    # move the solids, as they are created with centroids at (0,0,0)
    # 0.005 below due to the 0.01 increase in x-dir and z-dir grid lengths
    Fem.move(solid1,SVector{3,Float64}([ 30. + 0.005,  40. + 10., 30. + 0.005]))
    Fem.move(solid2,SVector{3,Float64}([ 30. + 0.005,  0.,        30. + 0.005]))
    #Fem.move(solid2,SVector{3,Float64}([ .5,  0., .5]))
    

    solids = [solid1 solid2]
    mats   = [material1 material2]

    #fixNodes(solid2,"bottom")
    #fixNodes(solid2,"left-right")


    c_dil     = sqrt(alumE/alumRho)
    dt        = grid.dy/c_dil
    dtime     = 0.1 * dt


    Tf      = 50e-3
    interval= 200
    #dtime   = 1.0e-8


    data                    = Dict()
    data["total_time"]      = Tf
    data["time"]            = 0.
    data["dt"]              = dtime
    data["dirichlet_grid"]  = [("back", (1,1,1)), 
                               ("front",(1,1,1)),
                               ("bottom",(1,1,1)),
                               ("left",(1,1,1)),
                               ("right",(1,1,1)),
                               ]      # => fix 

    data["friction"]        = 0.                             

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
    solve_explicit_dynamics_femp_3D_Contact(grid,solids,mats,basis,bodyforce,algo2,output2,fix,data)



    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))
# end

# @time main()
