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
using  WriteVTK

#function main()

    # problem parameters
	density       = 1050.
	E             = 1e6
	nu            = 0.3
	G             = 1.
	T             = 1.

	mu     = E/(2*(1+nu))
    lambda = E*nu/((1+nu)*(1-2*nu))

    Ro = 1.25
    Ri = 0.75
    l  = 2.8/2

    bspline = false

    xcenter   = 0.
    ycenter   = 0.

    # create the grid of a 1 x 1 square, with 20 x 20 cells
    if !bspline
       grid  =  Grid2D(-l, l, -l, l, 81, 81)
       basis = LinearBasis()
    else
       grid  =  Grid2D(0, 2*l, 0, 2*l, 101, 101)
       basis =  QuadBsplineBasis()
       xcenter   = l
       ycenter   = l
    end


    material  = NeoHookeanMaterial(E,nu,density)
    #material  = ElasticMaterial(E,nu,density,0,0)

	c_dil     = sqrt((material.lambda + 2*material.mu)/material.density)
	dt        = grid.dx/c_dil
	dtime     = 0.02 * dt


    #solid1   = FEM2D("vortex.msh")
    solid1   = FEM2D("vortex-regular.msh")
    

    # v0 = SVector{2,Float64}([vel  0.0])

    # # assign initial velocity for the particles
    # Fem.assign_velocity(solid1, v0)
    # Fem.assign_velocity(solid2,-v0)

    if (bspline) Fem.move(solid1,SVector{2,Float64}([ l,  l])) end
    
    fixNodes(solid1, "inner")
    fixNodes(solid1, "outer")

    solids = [solid1]
    mats   = [material]
    
 
    bodyforce = VortexBodyForce2D(G,T,Ri,Ro,density,mu,xcenter,ycenter)

    @printf("Total number of material points: %d \n", solid1.parCount)
    @printf("Total number of grid points:     %d\n", grid.nodeCount)
    @printf("Sound vel : %+.6e \n", c_dil)
    @printf("dt        : %+.6e \n", dtime)
    println(typeof(basis))
    println(solid1.vtk_cell)

    Tf      = 1.#5e-3 #3.5e-0
    interval= 50

	output2  = VTKOutput(interval,"vortex-femp-results/",["vx","sigmaxx"])
	fix      = EnergiesFix(solids,"vortex-femp-results/energies.txt")
    
    algo0    = USL(0.)
    algo1    = TLFEM(0.)
    algo2    = TLFEMFull(1e-8)
    
	plotGrid(output2,grid)
	plotParticles_2D(output2,solids,mats,0)
    #solve_explicit_dynamics_femp_2D(grid,solids,mats,basis,bodyforce,algo1,output2,fix,Tf,dtime)

 
    function dodo(solids)
           t = 0.
    counter = 0
    while t < Tf
        my_vtk_file = string("vortex-femp-results/","exact_","$(Int(counter))")
        vtmfile     = vtk_multiblock(my_vtk_file)

        g     = G * sin(π*t/T)

     
        solid = solids[1]
        xx    = solid.pos
        XX    = solid.pos0
        elems = solid.elems
        
        points = zeros(2,solid.nodeCount)       
        disp   = zeros(2,solid.nodeCount)
        
        cc = 1
        cells = MeshCell[]

        for ip=1:solid.nodeCount
            x     = XX[ip][1]-xcenter
            y     = XX[ip][2]-ycenter
            r     = sqrt(x*x+y*y)
            

            R     = (Ri+Ro)/2
            h     = 1 - 8*((r - R)/(Ri-Ro))^2 +16*((r - R)/(Ri-Ro))^4
            alpha = g*h

            points[1,ip] = (XX[ip][1]-xcenter)*cos(alpha) - (XX[ip][2]-ycenter)*sin(alpha) + xcenter
            points[2,ip] = (XX[ip][1]-xcenter)*sin(alpha) + (XX[ip][2]-ycenter)*cos(alpha) + ycenter
            
            disp[1,ip] = points[1,ip] - XX[ip][1]
            disp[2,ip] = points[2,ip] - XX[ip][2]
            #cc += 1
        end
        for e=1:solid.parCount
            inds =elems[e,:] 
            c    = MeshCell(VTKCellTypes.VTK_QUAD, inds)
            push!(cells, c)
        end

        vtkfile     = vtk_grid(vtmfile, points, cells)       
        vtkfile["Displacement", VTKPointData()] = disp
        outfiles    = vtk_save(vtmfile)

        t       += 0.001
        counter +=1
    end
end

#dodo(solids)


# end

# @time main()
