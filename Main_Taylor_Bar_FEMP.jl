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


# This file is for the Taylor anvil test presented in
# paper 'A generalized particle in cell metghod for explicit solid dynamics', V.P. Nguyen et all 2020.

push!(LOAD_PATH,"./")
#
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")


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

    A  = 0.065
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

    if vel == 190 
        # create the grid of a 1 x 1 square, with 20 x 20 cells
        L     = 16 
        grid  =  Grid3D(0, L, 0, 28, 0, L, 20, 60, 20) # for 190 m/s
    else
        L     = 40
        grid  =  Grid3D(0, L, 0, 28, 0, L, 50, 60, 50) # for 750 m/s
    end
    
   # basis = QuadBsplineBasis()
    basis = LinearBasis()

    # (E,nu,density,A,B,C,n,eps0dot,cellsize,parCount)
    
    #material  = NeoHookeanMaterial(E,nu,density)
    
    
	c_dil     = sqrt(E/density)
	dt        = grid.dy/c_dil
	dtime     = 0.1 * dt


    #solid1   = FEM3D("taylor-bar-47040.msh")
    #solid1   = FEM3D("taylor-bar-tet4.msh",material)
    solid1   = FEM3D("taylor-bar-52920.msh")

    
    material = JohnsonCookMaterial(E,nu,density,A,B,C,n,eps0dot,m,χ,Cp,.42,solid1.parCount)

    v0       = SVector{3,Float64}([0.0 -vel 0.])

    # assign initial velocity for the particles
    Fem.assign_velocity(solid1, v0)
    
    Fem.move(solid1,SVector{3,Float64}([ L/2, 1.6, L/2]))
    

    solids = [solid1]
    mats   = [material]
    

    @printf("Total number of material points: %d \n", solid1.parCount)
    @printf("Total mass: %+.6e \n", sum(solid1.mass))
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("Sound vel : %+.6e \n", c_dil)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))

    Tf      = 63e-3
    interval= 50

	output   = VTKOutput(interval,"taylor-femp-tet4-results/",["vx","sigmaxx"])
	fix      = EnergiesFix(solids,"taylor-femp-tet4-results/energies.txt")
    algo1    = USL(0.)
    algo2    = TLFEM(0.,1.)
    bodyforce = ConstantBodyForce3D([0., 0.,0.])


    data               = Dict()
    data["total_time"] = Tf
    data["dt"]         = dtime
    data["time"]       = 0.
    data["dirichlet_grid"] = [("bottom",(0,1,0))] # => fix bottom nodes on Y dir

	plotGrid(output,grid,0)
	plotParticles_3D(output,solids,mats,0)
    solve_explicit_dynamics_femp_3D(grid,solids,mats,basis,bodyforce,algo2,output,fix,data)

    D = abs(solid1.pos[1][1] - solid1.pos[3][1])
    L = abs(solid1.pos[2][2] - solid1.pos[6][2])
    #W = abs(solid1.pos[?][1] - solid1.pos[?][1]) # 5.08 from the bottom 
#  end

# @time main()
