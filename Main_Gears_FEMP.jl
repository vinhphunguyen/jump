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

# This file is for two gears problem.

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/juMP")
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
	g              = 0
	
	density = 7.85e-9;
	E       = 210e3;
	nu      = 0.3;
	fric    = 0.0;

	ω       = -4e3 # radian/sec (clockwise rotation)
	R       = 0.5;

    # create the grid of a 1 x 1 square, with 20 x 20 cells
	# and a basis: linear and CPDI-Q4 supported
    grid      =  Grid3D(-1.5,7,-3.5,3,0,0.25,65,65,4)
    basis     =  LinearBasis()


    solid1   = FEM3D("gear-A.inp")
    solid2   = FEM3D("gear-B.inp")

    material1 = ElasticMaterial(E,nu,density,solid1.parCount)
    material2 = ElasticMaterial(E,nu,density,solid2.parCount)

  
    
    # as the mesh was created with the center of the disk at (0,0)
	#move(solid1,SVector{2,Float64}([ 0.2+grid.dx  0.2+grid.dx]))
	#move(solid2,SVector{2,Float64}([ 0.8-grid.dx  0.8-grid.dx]))


	Fem.move(solid1,SVector{3,Float64}([0,0,.08]))
	Fem.move(solid2,SVector{3,Float64}([0.,0,.08]))


    solids = [solid1]#, solid2]
    mats   = [material1]#,material2]

    Tf       = 200e-6 #3.5e-0
    interval = 2
	dtime    = 0.2*grid.dx/sqrt(E/density)


	# data is a dictionary: 
	# data["pressure"]   = [(1,"tag1",f),(2,"tag2",g)]
	# data["total_time"] = Tf
	# data["dt"]         = dtime
	# data["time"]       = t

	function f(x,y,z,t)
		theta = atan(y,x) + ω * t 
		return (R*cos(theta),R*sin(theta),z,-R*ω*sin(theta),R*ω*cos(theta),-R*ω*sin(theta)*t,R*ω*cos(theta)*t)
	end

	data               = Dict()
	data["total_time"] = Tf
	data["dt"]         = dtime
	data["time"]       = 0.
    data["time_dirichlet_solid"] = [(1,"inner",f)]
    data["dirichlet_grid"] = [("front",(0,0,1)), 
                              ("back",(0,0,1))] # => fix bottom nodes on Y dir
                             

    # Symmetric BCs
	# fixYNodes(solid1, "fix")
	# fixYNodes(solid2, "fix")

	#output1  = PyPlotOutput(interval,"twodisks-results/","Two Disks Collision",(4., 4.))
	output2  = VTKOutput(interval,"gears-femp/",["pressure"])	
	fix      = EnergiesFix(solids,"gears-femp/energies.txt")

    algo1    = USL(0.)
    algo2    = TLFEM(0.,1.)

    body     = ConstantBodyForce3D(@SVector[0.,-g,0.])

	report(grid,solids,dtime)

    plotGrid(output2,grid,0)
    plotParticles_3D(output2,solids,mats,0)

	#reset_timer!
    solve_explicit_dynamics_femp_3D_Contact(grid,solids,mats,basis,body,fric,algo2,output2,fix,data)
    #print_timer()

	# #PyPlot.savefig("plot_2Disk_Julia.pdf")

# end

# @time main()
