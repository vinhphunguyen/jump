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


# Input file for the problem of a rubber block impacting a rigid wall
# Solved with the GPIC "Generalized in Cell" Method
# Output in folder "vertical-bar-drop-femp", with  vtk files and energies.txt

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

function main()

    # problem parameters
	g             = 0
	density       = 1050e-12
	youngModulus  = 1.0
	poissonRatio  = 0.3


	vel           = 50e3

    # create the grid of a 2x2x2 m3, with 15 x 15 x 15 cells
	# and a basis: linear 
    grid      =  Grid3D(0,2500,0,2500,0,2000,15,15,15)
    basis     =  LinearBasis()


    solid1    = FEM3D("bar8000.msh")
    material  = NeoHookeanMaterial(youngModulus,poissonRatio,density,solid1.parCount)

    
    # as the mesh was created with the center of the disk at (0,0)
	#move(solid1,SVector{2,Float64}([ 0.2+grid.dx  0.2+grid.dx]))
	#move(solid2,SVector{2,Float64}([ 0.8-grid.dx  0.8-grid.dx]))
	Fem.move(solid1,SVector{3,Float64}([2000/4,2000/4,2000/4]))

    Fem.assign_velocity(solid1, SVector{3,Float64}([0. -vel 0.0 ]))


    solids = [solid1]
    mats   = [material]

    Tf       = 0.08 #3.5e-0
    interval = 2
	dtime    = 0.1*grid.dx/sqrt(youngModulus/density)

	#output1  = PyPlotOutput(interval,"twodisks-results/","Two Disks Collision",(4., 4.))
	output2  = VTKOutput(interval,"vertical-bar-drop-femp/",["pressure"])
	fix      = DisplacementFemFix(solid1,"vertical-bar-drop-femp/",3)

    algo1    = USL(0.)
    algo2    = TLFEM(0.,1.)
    algo3    = TLFEMFull(0.)

    body     = ConstantBodyForce3D(@SVector[0.,-g,0.])


    data                    = Dict()
    data["total_time"]      = Tf
    data["dt"]              = dtime
    data["time"]            = 0.
    data["dirichlet_grid"]  = [("bottom",(0,1,0)),
                               ("left",(1,1,1)),
                               ("right",(1,1,1)),
                               ("front",(1,1,1)),
                               ("back",(1,1,1)),
                              ] # 
    #data["dirichlet_solid"] = [(1,"TopSurface",(0,1,0))] # => fix  nodes of 'TopSurface' group of solid 1 on Y dir
               

	report(grid,solids,dtime)

    plotGrid(output2,grid,0)
    #plotParticles_3D(output2,solids,0)

	#reset_timer!
    solve_explicit_dynamics_femp_3D(grid,solids,mats,basis,body,algo3,output2,fix,data)
    #print_timer()

	# #PyPlot.savefig("plot_2Disk_Julia.pdf")

end

@time main()
