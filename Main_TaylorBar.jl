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

# Input file for the Taylor bar problem: a metal bar impacting a rigid wall
# Solved with the standard material point method
# Output in folder "TaylorBar3D/", with vtu files and energies.txt

push!(LOAD_PATH,"./")
# import Gadfly
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")

# pyFig_RealTime = PyPlot.figure("MPM 2Disk Real-time",
#                                figsize=(8/2.54, 8/2.54), edgecolor="white", facecolor="white")

#include("./Grid.jl")
#include("./Problem.jl")

using Solid
using Grid
using Problem
using Output
using Algorithm
using Material
using BodyForce
using Basis
using Fix

#function main()

    # problem parameters
    E       = 115
    nu      = 0.31
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

    vel     = 190 # mm/s


	l0      = 25.4 #mm
	r0      = 3.8 #mm


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

    basis = LinearBasis()
    #basis = QuadBsplineBasis()

    ppc = 3
    fOffset = grid.dx/ppc

    coords1 = buildParticleForCylinder([L/2;  15; L/2], r0,
	                                   l0/2, fOffset, fOffset)

    
    material  = JohnsonCookMaterial(E,nu,density,A,B,C,n,eps0dot,m,χ,Cp,.42,length(coords1))
    solid1    = Solid3D(coords1,material)

    solid1.mass          .*= fOffset * fOffset * fOffset
    solid1.volume        .=  fOffset * fOffset * fOffset
    solid1.volumeInitial .=  fOffset * fOffset *fOffset

    v0 = SVector{3,Float64}([0.0 -vel 0])

    # assign initial velocity for the particles
    Solid.assign_velocity(solid1, v0)

    solids = [solid1]

    @printf("Total number of material points: %d \n", solid1.parCount)
    @printf("Total mass: %f \n", sum(solid1.mass))
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("Vol  : %+.6e \n", sum(solid1.volume))


    Tf      = 63e-3
    interval= 50

    
    c_dil     = sqrt(E/density)
    dt        = grid.dy/c_dil
    dtime     = 0.1 * dt

    data               = Dict()
    data["total_time"] = Tf
    data["dt"]         = dtime
    data["time"]       = 0.
    data["dirichlet_grid"] = [("bottom",(0,1,0))] # => fix bottom nodes on Y dir

    output2  = OvitoOutput(interval,"TaylorBar3D/",["pressure", "vonMises", "pstrain"])
    fix      = EmptyFix()
    algo1    = USL(0.)
    algo2    = MUSL(1.)

    plotGrid(output2,grid,0)

    plotParticles_3D(output2,solids,[grid.lx, grid.ly, grid.lz],
                        [grid.nodeCountX, grid.nodeCountY, grid.nodeCount],0)

    solve_explicit_dynamics_3D(grid,solids,basis,algo2,output2,fix,data)

# end

# @time main()
