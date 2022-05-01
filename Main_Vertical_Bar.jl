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
# Solved with the GPIC "Generalized Particle in Cell" Method
# Output in folder "vertical-bar/", with  vtk files and energies.txt

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


using Solid
using Grid
using Output
using Algorithm
using Material
using BodyForce
using Basis
using Fix
using Util
using Problem

function main()

   # problem parameters
	g             = 0#1e6
	density       = 1050e-12
	youngModulus  = 1.0
	poissonRatio  = 0.3

	vel           = 50e3

   # create the grid of a 1 x 1 square, with 20 x 20 cells
	# and a basis: linear 
   grid      =  Grid3D(0,2000,0,2000,0,2000,15,15,15)
   basis     =  LinearBasis()

   ppc      = 2
   fOffset  = dx = grid.dx/ppc
   coords   = buildParticleForBlock([1000.;  1000; 1000.], 1000.,1000.,1000., fOffset)
   material = NeoHookeanMaterial(youngModulus,poissonRatio,density,length(coords))

   solid1   = Solid3D(coords,material)

   solid1.mass          .*= dx * dx
   solid1.volume        .= dx * dx
   solid1.volumeInitial .= dx * dx
    

   Solid.assign_velocity(solid1, SVector{3,Float64}([0. -vel 0.0 ]))


   solids = [solid1]

   Tf       = 0.1 #3.5e-0
   interval = 2
	dtime    = 0.08*grid.dx/sqrt(youngModulus/density)

	#output1  = PyPlotOutput(interval,"twodisks-results/","Two Disks Collision",(4., 4.))
	output2  = OvitoOutput(interval,"vertical-bar/",["pressure"])
	fix      = EmptyFix()#
	fix      = EmptyFix()

   algo1    = USL(0.)
   algo2    = MUSL(0.99)

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
                              ] # => fix bottom nodes on Y dir

	report(grid,solids,dtime)

   plotGrid(output2,grid,0)
   plotParticles_3D(output2,solids,[grid.lx, grid.ly, grid.lz],
		            [grid.nodeCountX, grid.nodeCountY, grid.nodeCount],0)

	#reset_timer!
   solve_explicit_dynamics_3D(grid,solids,basis,algo2,output2,fix,data)
    
    #print_timer()

	# #PyPlot.savefig("plot_2Disk_Julia.pdf")

end

@time main()
