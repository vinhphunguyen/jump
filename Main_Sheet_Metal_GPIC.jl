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

using DelimitedFiles
using Velocity

function main()

    # elastic parameters
    alumRho   = 2.71e-6
    alumE     = 69
    alumNu    = 0.33

    # Johnson-Cook parameters
    A       = 0.379
    B       = 0.163
    C       = 13.05
    Tm      = 1600.
    χ       = 0.
    Cp      = 0.

    fric    = 0.0

    # grid creation: 200x10x200 mm3 with origin at 0,0,0 by 100x20x100 cells
    # 10: thickness dir, i do not know how much the sheet will deform!!!
    # the actual grid length in x-dir is 2.01 so that there is no FE nodes exactly locating on the grid line!!!
    # similar for z-dir
    # See Basis.jl, getNearestGridPoints(points,start,xp,grid::Grid3D), floor(Int64,(xp[1]-grid.xmin) * fLength_Cell_x + 1.)
    # linear basis for the grid
    grid  =  Grid3D(0.,200.01, -10.,15., 0, 200.01, 101, 31, 101)
    basis = LinearBasis()

    solid1   = FEM3D("sheet-metal-indenter.msh")
    solid2   = FEM3D("sheet-metal-blanket.msh")
    

    material1 = RigidMaterial(alumRho)
    material2 = VoceMaterial(alumE,alumNu,alumRho,A,B,C,χ,Cp,.42,solid2.parCount) 

    Fem.move(solid1,SVector{3,Float64}([116.953+0.005,12-1.35-2.8+0.0011,12.911+0.005])) # 100 in the x to move the indenter to the correct pos: do not change
                                                              # 12.7+1.5: to move the indenter in the y dir: 1.5 = thickness of the sample, 12.7 is the indenter's diameter
                                                              # 100+0.005 (0.005=0.01/2 as the z length is 0.51 slighyl larger than 0.5) is to move it to the middle (z-dir): do not change
    Fem.move(solid2,SVector{3,Float64}([ 0.005,  0., 0.005])) # 5 in y dir, to move all objects up to 5 (will change it according to the deformation)
    #Fem.move(solid2,SVector{3,Float64}([ .5,  0., .5]))
    

    solids = [solid1, solid2]
    mats   = [material1, material2]


    c_dil     = sqrt(alumE/alumRho)
    dt        = grid.dy/c_dil
    dtime     = 0.09 * dt
    dt_factor = 0.5


    Tf      = 104.0
    interval= 20

    vT = readdlm("velocity_ISF.txt", '\t', Float64, '\n', skipstart=1)
    velo_data = VelocityData(vT[:,1],vT[:,2],vT[:,3],vT[:,4])

    data                    = Dict()
    data["total_time"]      = Tf
    data["time"]            = 0.
    data["dt"]              = dtime
    data["dt_factor"]       = dt_factor
    data["friction"]        = fric
    data["dirichlet_grid"]  = [("back", (1,1,1)), 
                               ("front",(1,1,1)),
                               ("left",(1,1,1)),
                               ("right",(1,1,1)),
                               ]      # => fix bottom nodes on X/Y/Z     dir             

    #data["dirichlet_solid"] = [(2,"symmetry",(1,1,1))] # => fix  nodes of 'fix' group of solid 2 on Y dir
#                                                                            
    data["rigid_body_velo_data"] = Array{Tuple{Int64,VelocityData},1}([(1,velo_data)])            # => solid 1 has a velo given by velo_func


    output2  = VTKOutput(interval,"sheet-metal-results/",["vx","sigmaxx"])
    fix      = EmptyFix() #(solid2,"bottom","scratch-3D-results/forces.txt")
    algo1    = USL(0.)
    algo2    = TLFEM(0.,0.999)
    bodyforce = ConstantBodyForce3D([0.,0.,0.])



    @printf("Total number of material points: %d \n", solid1.parCount)#+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))

    plotGrid(output2,grid,0)
    plotParticles_3D(output2,solids,mats,0)

    #reset_timer!()
    #solve_explicit_dynamics_femp_3D(grid,solids,mats,basis,bodyforce,algo2,output2,fix,data)
    #print_timer()
    solve_explicit_dynamics_femp_3D_Contact(grid,solids,mats,basis,bodyforce,algo2,output2,fix,data)


    @printf("Total number of material points: %d \n", solid1.parCount)#+solid2.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))
end

@time main()
