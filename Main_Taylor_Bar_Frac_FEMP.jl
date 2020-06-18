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
	

    E   = 200
    nu  = 0.33
    density = 7.85e-06

    A  = 490e-3
    B       = 0.383
    C       = 0.0123
    n       = 0.45
    eps0dot = 5e-4
    Tr      = 293
    Tm      = 1800
    m       = 0.94
    χ       = 0.9
    Cp      = 452

    D1      = 0.0705
    D2      = 1.732
    D3      = -0.54
    D4      = 0.015
    D5      = 0.0

	vel     = 600# mm/s

    L     = 40
    grid  =  Grid3D(0, L, 0, 32, 0, L, 50, 60, 50) # for 750 m/s
    
    
   # basis = QuadBsplineBasis()
    basis = LinearBasis()

    # (E,nu,density,A,B,C,n,eps0dot,cellsize,parCount)
    
    #material  = NeoHookeanMaterial(E,nu,density)
    
    
	c_dil     = sqrt(E/density)
	dt        = grid.dy/c_dil
	dtime     = 0.1 * dt


    #solid1   = FEM3D("taylor-bar-47040.msh")
    #solid1   = FEM3D("taylor-bar-tet4.msh",material)
    solid1   = FEM3D("taylor-bar-frac.inp")

    strength  = JohnsonCookMaterial(E,nu,density,A,B,C,n,eps0dot,m,χ,Cp,.42,solid1.parCount)
    damage    = JohnsonCookDamage(D1,D2,D3,D4,D5,eps0dot,Tr,Tm,solid1.parCount)
    material  = JohnsonCookMaterialWithDamage(strength,damage)

    v0       = SVector{3,Float64}([0.0 -vel 0.])

    # assign initial velocity for the particles
    Fem.assign_velocity(solid1, v0)
    
    Fem.move(solid1,SVector{3,Float64}([ L/2, 1.6, L/2]))
    

    solids = [solid1]
    mats   = [material]

   
 
    fixYForBottom(grid)
    

    @printf("Total number of material points: %d \n", solid1.parCount)
    @printf("Total mass: %+.6e \n", sum(solid1.mass))
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("Sound vel : %+.6e \n", c_dil)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))

    Tf      = 63e-3
    interval= 50

	output   = VTKOutput(interval,"taylor-femp-frac-results/",["vx","sigmaxx"])
	fix      = EnergiesFix(solids,"taylor-femp-frac-results/energies.txt")
    algo1    = USL(0.)
    algo2    = TLFEM(0.,1.)
    bodyforce = ConstantBodyForce3D([0., 0.,0.])

	plotGrid(output,grid,0)
	plotParticles_3D(output,solids,mats,0)
    solve_explicit_dynamics_femp_3D(grid,solids,mats,basis,bodyforce,algo2,output,fix,Tf,dtime)

    #W = abs(solid1.pos[?][1] - solid1.pos[?][1]) # 5.08 from the bottom 
#  end

# @time main()
