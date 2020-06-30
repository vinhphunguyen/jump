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

    # elastic parameters
    alumRho   = 8.94e-6
    alumE     = 115
    alumNu    = 0.31

    # Johnson-Cook parameters
    A       = 0.145
    B       = 0.0
    C       = 0#0.013
    n       = 0.
    eps0dot = 1e-3
    Tm      = 1600.
    m       = 0.
    χ       = 0.
    Cp      = 0.

    fric    = 0.

	# grid creation
	grid  = Grid2D(0.,2.01, 0.,1., 51, 61)
	basis = LinearBasis()

    solid1   = FEM3D("SPHERE1.inp")
    solid2   = FEM3D("SUBSTRATE1.inp")
    

    material1 = RigidMaterial(alumRho)
    material2 = JohnsonCookMaterial(alumE,alumNu,alumRho,A,B,C,n,eps0dot,m,χ,Cp,.42,solid2.parCount) 
    #ElastoPlasticMaterial(alumE,alumNu,alumRho,alumFy,alumK,solid2.parCount)

    Fem.move(solid1,SVector{3,Float64}([ 0.005,0.05,0.025]))
    Fem.move(solid2,SVector{3,Float64}([ 0.005,  0., 0.025]))
    #Fem.move(solid2,SVector{3,Float64}([ .5,  0., .5]))
    

    solids = [solid1, solid2]
    mats   = [material1, material2]


    c_dil     = sqrt(alumE/alumRho)
    dt        = grid.dy/c_dil
    dtime     = 0.3 * dt


    Tf      = 12e-3
    interval= 20
    
    function velo_func(t)
        if t < 1.5e-3
            vx = 0.
            vy = -50.            
        else
            vx = 50.
            vy = 0.            
        end
        return (vx,vy)
    end

    data                    = Dict()
    data["total_time"]      = Tf
    data["time"]            = 0.
    data["dt"]              = dtime
    data["friction"]        = fric
   

    data["dirichlet_solid"] = [(2,"symmetry",(0,0,1))] # => fix  nodes of 'fix' group of solid 2 on Y dir
#                                                                            
    data["rigid_body_velo"] = [(1,velo_func)]            # => solid 1 has a velo given by velo_func                             
                             


	output2  = VTKOutput(interval,"scratch-2D-results/",["vx","sigmaxx"])
	fix      = ReactionForceFix(solid2,"scratch-2D-results/forces.txt")
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
