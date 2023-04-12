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

# Input file for the vibration of a 3D compliant cantilever beam 
# Solved with the GPIC method, 
# Output in folder "cantilever-box-femp/", with lammps dump files and energies.txt

push!(LOAD_PATH,"./")

# import Gadfly
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")
using TimerOutputs
# pyFig_RealTime = PyPlot.figure("MPM 2Disk Real-time",
#                                figsize=(8/2.54, 8/2.54), edgecolor="white", facecolor="white")

#include("./Grid.jl")
#include("./Problem.jl")

using Fem
using Solid
using Grid
using FemMPM
using Output
using Algorithm
using Material
using BodyForce
using Basis
using Fix
using Util

#function main()

    # problem parameters
	g             = 10.
	density       = 1050
	youngModulus  = 1e6
	poissonRatio  = 0.3

    # create the grid of a 8 x 8 x 2, with 24 x 24 x 2 cells
	# and a basis: 
    grid      =  Grid3D(0,8,0,8,0,2,25,25,7)
    basis     =  LinearBasis()


    #solid1   = FEM3D("box-beam.msh") # 8-node hexahedron 
    solid1   = FEM3D("box-beam-tet4.msh")
    material = NeoHookeanMaterial(youngModulus,poissonRatio,density,solid1.parCount)

    
    # as the mesh was created with the center of the disk at (0,0)
	#move(solid1,SVector{2,Float64}([ 0.2+grid.dx  0.2+grid.dx]))
	#move(solid2,SVector{2,Float64}([ 0.8-grid.dx  0.8-grid.dx]))
	Fem.move(solid1,SVector{3,Float64}([0.0,6.,0.5]))

    #Fem.assign_velocity(solid1, SVector{3,Float64}([0. -1000. 0.0 ]))

    solids = [solid1]
    mats   = [material]

    Tf       = 3.0 #3.5e-0
    interval = 20
	dtime    = 0.1*grid.dx/sqrt(youngModulus/density)
    dtime    = 1e-5 #0.2*0.08/sqrt(youngModulus/density)

	#output1  = PyPlotOutput(interval,"twodisks-results/","Two Disks Collision",(4., 4.))
	output2  = VTKOutput(interval,"cantilever-box-femp/",["pressure","sigmaxx"])
	fix      = DisplacementFemFix(solid1,"cantilever-box-femp/",2)

    algo1    = USL(0.)
    algo2    = TLFEM(0.,1.)

    body     = ConstantBodyForce3D(@SVector[0.,-g,0.])

    data                    = Dict()
    data["total_time"]      = Tf
    data["dt"]              = dtime
    data["time"]            = 0.
    data["dirichlet_solid"]  = [(1,"fix",(1,1,1))] 

	report(grid,solids,dtime)

    plotGrid(output2,grid,0)
    plotParticles_3D(output2,solids,mats,0)

	#reset_timer!
    solve_explicit_dynamics_femp_3D(grid,solids,mats,basis,body,algo2,output2,fix,data)
    #print_timer()

	# #PyPlot.savefig("plot_2Disk_Julia.pdf")
# end

# @time main()
