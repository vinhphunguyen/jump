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


# Taylor anvil test with axis-symmetric elements

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
	

    E   = 115
    nu  = 0.31
    density = 8.94e-06

    A       = 0.065
    B       = 0.356
    C       = 0#0.013
    n       = 0.37
    eps0dot = 1e-3
    Tm      = 1600
    m       = 0.
    χ       = 0.
    Cp      = 0.

	vel     = 190# mm/s
    #vel     = 750# mm/s

 
    L     = 8 
    grid  =  Grid2D(0, L, 0, 28, 20, 60) # for 190 m/s
  
    
   # basis = QuadBsplineBasis()
    basis = LinearBasis()

    # (E,nu,density,A,B,C,n,eps0dot,cellsize,parCount)
    
    #material  = NeoHookeanMaterial(E,nu,density)
    
    
	c_dil     = sqrt(E/density)
	dt        = grid.dy/c_dil
	dtime     = 0.1 * dt


    #solid1   = FEM3D("taylor-bar-47040.msh")
    #solid1   = FEM3D("taylor-bar-tet4.msh",material)
    solid1   = FEMAxis("taylor-bar-axis.msh")


    material = JohnsonCookMaterial(E,nu,density,A,B,C,n,eps0dot,m,χ,Cp,.42,solid1.parCount)

    v0       = SVector{2,Float64}([0.0 -vel])

    # assign initial velocity for the particles
    Fem.assign_velocity(solid1, v0)
    
    Fem.move(solid1,SVector{2,Float64}([ 0.05, 1.6]))
    

    solids = [solid1]
    mats   = [material]

   
 
    fixXForLeft(grid)
    #fixXForLeft(solid1)
    fixYForBottom(grid)
    

    @printf("Total number of material points: %d \n", solid1.parCount)
    @printf("Total mass: %+.6e \n", sum(solid1.mass))
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("Sound vel : %+.6e \n", c_dil)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))

    Tf      = 63e-3
    interval= 50

	output   = VTKOutput(interval,"taylor-femp-axis-results/",["vx","sigmaxx"])
	fix      = EnergiesFix(solids,"taylor-femp-axis-results/energies.txt")
    algo1    = USL(0.)
    algo2    = TLFEM(0.,1.)
    bodyforce = ConstantBodyForce2D([0., 0.])

	plotGrid(output,grid)
	plotParticles_2D(output,solids,mats,0)
    solve_explicit_dynamics_femp_3D(grid,solids,mats,basis,bodyforce,algo2,output,fix,Tf,dtime)

    D = 2*abs(solid1.pos[1][1] - solid1.pos[2][1])
    L =   abs(solid1.pos[1][2] - solid1.pos[4][2])
    @printf("Numerical diameter: %f \n", D)
    @printf("Numerical length: %f \n", L)
    #W = abs(solid1.pos[?][1] - solid1.pos[?][1]) # 5.08 from the bottom 
#  end

# @time main()
