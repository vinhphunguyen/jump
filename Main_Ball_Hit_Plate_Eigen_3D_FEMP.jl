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

    # problem parameters [mm,ms,kN,kg, GPa]
  
    v         = 500.0

    steelRho   = 2200e-9
    steelE     = 115. # GPa
    steelNu    = 0.3

    alumRho   = 2200e-9
    alumE     = 94.5 # GPa
    alumNu    = 0.3
    Gf1        = 2665e-4 # J/m2 N/m = 10-3 kN / 10^3 mm
    Gf2        = 4.3e-6 # J/m2 N/m = 10-3 kN / 10^3 mm



	# grid creation
	grid  = Grid3D(-40,40.0, 0.,24.0, -40,40, 71, 11, 71)
	basis = LinearBasis()

    solid1   = FEM3D("sphere-eigen.msh")
    solid2   = FEM3D("cylinder-eigen.msh")
    


    material1 = ElasticMaterial(steelE,steelNu,steelRho,Gf1,0.1,solid1.parCount)
    material2 = ElasticMaterial(alumE,alumNu,alumRho,Gf2,0.1,solid2.parCount) 
    #ElastoPlasticMaterial(alumE,alumNu,alumRho,alumFy,alumK,solid2.parCount)

    v0 = SVector{3,Float64}([0.0 -v 0.])

    # assign initial velocity for the particles
    Fem.assign_velocity(solid1, v0)
    Fem.move(solid1,SVector{3,Float64}([ 0., 17.5, 0.]))
    Fem.move(solid2,SVector{3,Float64}([ 0.,  2., 0.]))
    #Fem.move(solid2,SVector{3,Float64}([ .5,  0., .5]))
    

    solids = [solid1 solid2]
    mats   = [material1 material2]


    c_dil     = sqrt(steelE/steelRho)
    dt        = grid.dx/c_dil
    dtime     = 0.56 * dt


    Tf      = 100e-3
    interval= 40
    #dtime   = 1.0e-8

    data                    = Dict()
    data["total_time"]      = Tf
    data["time"]            = 0.
    data["dt"]              = dtime
    data["alpha"]           = 1.0
    data["kappa"]           = 0.001 # 1.25h where h is element size
    # data["friction"]        = fric
    # data["dirichlet_grid"]  = [("back", (0,0,1)), 
    #                            ("front",(1,1,1)),
    #                            ("bottom",(1,1,1)),
    #                            ("left",(1,1,1)),
    #                            ("right",(1,1,1)),
    #                            ]      # => fix bottom nodes on X/Y/Z     dir             

    data["dirichlet_solid"] = [(2,"fix",(1,1,1))] # => fix  nodes of 'fix' group of solid 2 on all dir
#                                                                            
    #data["rigid_body_velo"] = [(1,velo_func)]            # => solid 1 has a velo given by velo_func       


	output2  = VTKOutput(interval,"eigen-ball-plate-results/",["vx","sigmaxx"])
	fix      = EnergiesFix(solids,"eigen-ball-plate-results/energies.txt")
    algo1    = USL(0.)
    algo2    = TLFEM(0.,0.999)
    bodyforce = ConstantBodyForce3D([0.,0.,0.])



    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))

	plotGrid(output2,grid,0)
	plotParticles_3D(output2,solids,mats,0)
    solve_explicit_dynamics_femp_eigen_erosion_3D(grid,solids,mats,basis,bodyforce,algo2,output2,fix,data)


    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))
# end

# @time main()
