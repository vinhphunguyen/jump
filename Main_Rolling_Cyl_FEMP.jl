
# Phu Nguyen, Monash University
# 20 March, 2020 (Coronavirus outbreak)

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
using Mesh

#function main()

    # problem parameters
	# sphere
	K1      = 7;            # bulk modulus
	mu1     = 1.5;          # shear modulus
	lambda1 = K1 - 2/3*mu1;
	E1      = 9*K1*mu1/(3*K1+mu1);
	rho1    = 1e-9;           # density
	nu1     = (3*K1-2*mu1)/2/(3*K1+mu1);

	# plane = 10 times harder
	K2      = 70;          # bulk modulus
	mu2     = 15;          # shear modulus
	lambda2 = K2 - 2/3*mu2;
	E2      = 9*K2*mu2/(3*K2+mu2);
	rho2    = 1e-8;         # density
	nu2     = (3*K2-2*mu2)/2/(3*K2+mu2);

	theta = pi/3; # inclination angle
	fric  = 0.6;  # friction coefficient
	g     = 10e3; # gravity
	bodyf = [g*sin(theta) -g*cos(theta)];

    # create the grid of a 1 x 1 square, with 20 x 20 cells
	# and a basis: linear and CPDI-Q4 supported
	ra = 0.5e3;
	l  = 1.9e3;
	h  = 0.25e3;
	w  = h + 2*ra;
	ratio = l/w;

	noX0      = 50;        # number of elements along X direction
	noY0      = noX0/ratio; # number of elements along Y direction (to have square cells)
    grid      =  Grid2D(0,l, 0, w, noX0+1, floor(Int64,noY0+1))
    basis     =  LinearBasis()

 
    #material = NeoHookeanMaterial(youngModulus,poissonRatio,density)

    solid1   = FEM2D("rolling-cyl-circle.msh")
    solid2   = Rigid2D("rolling-cyl-plane.msh")

    material1 = ElasticMaterial(E1,nu1,rho1,solid1.parCount)
    material2 = RigidMaterial(rho2)
    #material2 = ElasticMaterial(E2,nu2,rho2,solid2.parCount)
    

    # as the mesh was created with the center of the disk at (0,0)
	#move(solid1,SVector{2,Float64}([ 0.2+grid.dx  0.2+grid.dx]))
	#move(solid2,SVector{2,Float64}([ 0.8-grid.dx  0.8-grid.dx]))
	Fem.move(solid1,SVector{2,Float64}([ ra+50,  ra+110]))
	Fem.move(solid2,SVector{2,Float64}([ 10.,  2.]))

	# initial stress due to gravity
	# funcs = zeros(4)
	# ders  = zeros(2,4)
	# for p = 1:solid1.parCount
	# 	xx = solid1.pos
	#     elemNodes  =  @view solid1.elems[p,:]  
	# 	elemNodes0 =        solid1.elems[p,:]  
	# 	coords     =  @view xx[elemNodes]
    
	# 	J          = lagrange_basis_derivatives!(funcs, ders, Quad4(), [0.,0.], coords)

 #        # only Quad4, but we do not want to use Tri3, do we?
	# 	yp         = funcs[1] * coords[1][2] + funcs[2] * coords[2][2] + funcs[3] * coords[3][2] + funcs[4] * coords[4][2]

	#     hh                = 130. +2*ra-yp;
	#     solid1.stress[p][1,1] = -rho1*g*hh;
	#     solid1.stress[p][2,2] = -rho1*g*hh;
	# end

    solids = [solid1, solid2]
    mats   = [material1,material2]

    Tf       = 0.3
    interval = 100
	
	c_dil     = sqrt(E2/rho2)
	dt        = grid.dy/c_dil
	dtime     = 0.1 * dt


	#output1  = PyPlotOutput(interval,"twodisks-results/","Two Disks Collision",(4., 4.))
	output2  = VTKOutput(interval,"rolling-cyl-femp/",["pressure"])
	fix      = CenterOfMassFix(solid1,"rolling-cyl-femp/numerical.txt")

    algo2    = USL(1e-9)
    algo1    = TLFEM(1e-9,1.0)
    body     = ConstantBodyForce2D(bodyf)

	report(grid,solids,dtime)

    plotGrid(output2,grid)
    plotParticles_2D(output2,solids,mats,0)

    # fixXForBottom(grid)
    # fixYForBottom(grid)
    # fixYNodes(solid2, "BotSurface")

    data                    = Dict()
    data["total_time"]      = Tf
    data["time"]            = 0.
    data["dt"]              = dtime
    data["friction"]        = fric
    #data["dirichlet_grid"]  = [("bottom",(1,1)),]      # => fix bottom nodes on X/Y/Z     dir   
    data["rigid_body_velo"] = [(2,t -> (0.,0.))]            # => solid 1 has a velo given by velo_func            

    #data["dirichlet_solid"] = [(2,"symmetry",(0,0,1))] # => fix  nodes of 'fix' group of solid 2 on Y dir
#                                                                           

	#reset_timer!
    solve_explicit_dynamics_femp_2D_Contact(grid,solids,mats,basis,body,algo1,output2,fix,data)
    #print_timer()
	# plotting energies
	# analytical solution
	# filename = string("rolling-cyl-femp/","exact-stick.txt")
 #    if (isfile(filename)) rm(filename) end
 #    file = open(filename, "a")
 #    tt = LinRange(0.,Tf,30);
	
	# for i=1:length(tt)
	# 	#xx1 = ra + 0.5*g*tt[i]^2*(sin(theta)-fric*cos(theta));
	# 	xx1 = ra + 1/3*g*tt[i]^2*sin(theta);
	#     write(file, "$(tt[i]) $(xx1)\n")
 #    end
 #    close(file)

# end

# @time main()
