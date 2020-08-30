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

push!(LOAD_PATH,"./")
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

# elastic parameters
alumRho   = 8.94e-6
alumE     = 115
alumNu    = 0.31

# Johnson-Cook parameters
A       = 0.145
B       = 0.23
C       = 0.013
n       = 0.
eps0dot = 1e-3
Tm      = 1600.
m       = 0.
χ       = 0.
Cp      = 0.

fric    = 0

# grid creation: 2x1x0.5 with origin at 0,0,0 by 80x50x20 cells
# the actual grid length in x-dir is 2.01 so that there is no FE nodes exactly locating on the grid line!!!
# similar for z-dir
# See Basis.jl, getNearestGridPoints(points,start,xp,grid::Grid3D), floor(Int64,(xp[1]-grid.xmin) * fLength_Cell_x + 1.)
# linear basis for the grid
grid  =  Grid3D(0.,2.01, 0.,1., 0, 0.51, 81, 51, 21)
#grid  = Grid3D(0.,2.1, 0.,1., -0.025, 0.525, 81, 51, 23)
basis = LinearBasis()
solid1   = FEM3D("sphere.msh")
solid2   = FEM3D("SUBSTRATE1.inp")
#solid2   = FEM3D(“SUBSTRATE2.inp”) # non-conforming mesh
material1 = RigidMaterial(alumRho)
material2 = JohnsonCookMaterial(alumE,alumNu,alumRho,A,B,C,n,eps0dot,m,χ,Cp,.42,solid2.parCount)
#ElastoPlasticMaterial(alumE,alumNu,alumRho,alumFy,alumK,solid2.parCount)
Fem.move(solid1,SVector{3,Float64}([ 0.5+0.005,0.15+0.5+1.2*grid.dy,0.25+0.005])) # 0.5 in the x to move the indenter to the correct pos: do not change
# 0.5+0.15+ to move the indenter in the y dir: 0,5 = height of the sample, 0.15 is the indenter’s radius
# 0.25+0.005 (0.005=0.01/2 as the z length is 0.51 slighyl larger than 0.5) is to move it to the middle (z-dir): do not change
Fem.move(solid2,SVector{3,Float64}([ 0.005,  0., 0.005]))
#Fem.move(solid2,SVector{3,Float64}([ .5,  0., .5]))    

solids = [solid1, solid2]
mats   = [material1, material2]


c_dil     = sqrt(alumE/alumRho)
dt        = grid.dy/c_dil
dtime     = 0.09 * dt
dt_factor = 0.5

Tf      = 12e-3
interval= 50

function velo_func(t)
    if t < 0.5e-3
        vx = 0.
        vy = -50.
        vz = 0.
    else
        vx = 50.
        vy = 0.
        vz = 0.
    end
    return (vx,vy,vz)
end

data                    = Dict()
data["total_time"]      = Tf
data["time"]            = 0.
data["dt"]              = dtime
data["dt_factor"]       = dt_factor
data["friction"]        = fric
data["dirichlet_grid"]  = [("back", (1,1,1)), 
                           ("front",(1,1,1)),
                           ("bottom",(1,1,1)),
                           ("left",(1,1,1)),
                           ("right",(1,1,1)),
                           ]      # => fix bottom nodes on X/Y/Z     dir             

#data["dirichlet_solid"] = [(2,"symmetry",(1,1,1))] # => fix  nodes of 'fix' group of solid 2 on Y dir
#                                                                            
data["rigid_body_velo"] = Array{Tuple{Int64,Function},1}([(1,velo_func)])            # => solid 1 has a velo given by velo_func                             



output2  = VTKOutput(interval,"scratch-3D-results/",["vx","sigmaxx"])
fix      = ReactionForceFix(solid2,"bottom","scratch-3D-results/forces.txt")
algo1    = USL(0.)
algo2    = TLFEM(0.,0.999)
bodyforce = ConstantBodyForce3D([0.,0.,0.])



@printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
@printf("Total number of grid points:     %d\n", grid.nodeCount)
@printf("dt        : %+.6e \n", dtime)
println(typeof(basis))

plotGrid(output2,grid,0)
plotParticles_3D(output2,solids,mats,0)

#reset_timer!()
solve_explicit_dynamics_femp_3D_Contact(grid,solids,mats,basis,bodyforce,algo2,output2,fix,data)
#print_timer()
#solve_explicit_dynamics_femp_3D_Contact(grid,solids,mats,basis,bodyforce,algo2,output2,fix,data)


@printf("Total number of material points: %d \n", solid1.parCount+solid2.parCount)
@printf("Total number of grid points:     %d\n", grid.nodeCount)
@printf("dt        : %+.6e \n", dtime)
println(typeof(basis))
# end

# @time main()
