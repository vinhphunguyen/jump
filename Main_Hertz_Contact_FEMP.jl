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

# Hertzian contact problem:
# two cylinders in contact, 3D
# taken from:Hertmann, CMAME, contact domain method part II

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/juMP")
#
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")
using TimerOutputs
using FunctionWrappers

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
    # units: mm,s,N,MPa
    # elastic parameters
    alumRho   = 8.94e-9
    alumE     = 200.
    alumNu    = 0.3


    p0      = 0.0625 # N/mm
    fric    = 0.0

	# grid creation
	grid  = Grid3D(-8.4,8.4, -8.0,8.4, 0, 1.1, 72, 121, 2)
	basis = LinearBasis()

    # solid1   = FEM3D("disk-half-bot.msh")
    # solid2   = FEM3D("disk-half-up.msh")
    # Fem.move(solid1,SVector{3,Float64}([ 0.2,  -8, 0.05]))
    # Fem.move(solid2,SVector{3,Float64}([ 0.2,  8.0, 0.05]))


    solid1   = FEM3D("disk-half-bot.inp")
    solid2   = FEM3D("disk-half-up.inp")
    

    material1 = ElasticMaterial(alumE,alumNu,alumRho,solid1.parCount) 
    material2 = ElasticMaterial(alumE,alumNu,alumRho,solid2.parCount) 
    #ElastoPlasticMaterial(alumE,alumNu,alumRho,alumFy,alumK,solid2.parCount)

    #Fem.move(solid1,SVector{3,Float64}([ 0.5,0.16+0.5,0.25]))
 

    #Fem.assign_velocity(solid2,[0,-500,0])
    

    solids = [solid1, solid2]
    mats   = [material1, material2]


    c_dil     = sqrt(alumE/alumRho)
    dt        = grid.dy/c_dil
    dtime     = 0.2 * dt


    Tf      = dtime
    t0      = 0.0001  #ms
    interval= 20
    
    function force_func(t)
        if t < t0
            px = (p0/t0)*t
        else
            px = p0
        end
        return (px)
    end

    data                    = Dict()
    data["total_time"]      = Tf
    data["time"]            = 0.
    data["dt"]              = dtime
    data["friction"]        = fric
    data["dirichlet_grid"]  = [
                               ("bottom",(1,1,1)),
                               ("front",(0,0,1)),
                               ("back",(0,0,1))
                               ]      # => fix bottom nodes on X/Y/Z     dir             

    data["dirichlet_solid"] = [(1,"bottom",(1,1,1))] # => fix  nodes of 'fix' group of solid 2 on Y dir
    data["pressure"]   = Array{Tuple{Int64,String,Function},1}([(2,"force",force_func)])
#                                                                                                     
                             

	output2  = VTKOutput(interval,"Hertz-results/",["vx","sigmaxx"])
	fix      = EmptyFix()#ReactionForceFix(solid1,"bottom", "Hertz-results/forces.txt")
    algo2    = TLFEM(0.,0.999)
    bodyforce = ConstantBodyForce3D([0.,0.,0.])



    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))

	plotGrid(output2,grid,0)
	plotParticles_3D(output2,solids,mats,0)

    reset_timer!()
    solve_explicit_dynamics_femp_3D_Contact(grid,solids,mats,basis,bodyforce,algo2,output2,fix,data)
    print_timer()
    #solve_explicit_dynamics_femp_3D_Contact(grid,solids,mats,basis,bodyforce,algo2,output2,fix,data)

    
    fix1 = SpatialReactionForceFix(solid1,"boundary", "Hertz-results/forces.txt")
    compute_femp(fix1,0.)
    closeFile(fix1)

    @printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))
# end

# @time main()
