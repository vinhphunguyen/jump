
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
	fric  = 0.0;  # friction coefficient
	g     = 10e3; # gravity
	bodyf = [g*sin(theta) -g*cos(theta)];
	bodyf = [0. 0.];

    # create the grid of a 1 x 1 square, with 20 x 20 cells
	# and a basis: linear and CPDI-Q4 supported
	ra = 0.5e3;
	l  = 1.65e3;
	h  = 0.25e3;
	w  = h + 2*ra;
	ratio = l/w;

	noX0      = 30;        # number of elements along X direction
	noY0      = noX0/ratio; # number of elements along Y direction (to have square cells)
    grid      =  Grid2D(0,l, 0, w, noX0+1, floor(Int64,noY0+1))
    basis     =  LinearBasis()

    material1 = ElasticMaterial(E1,nu1,rho1,0,0)
    material2 = ElasticMaterial(E2,nu2,rho2,0,0)
    #material = NeoHookeanMaterial(youngModulus,poissonRatio,density)

    solid1   = FEM2D("rolling-cyl-circle.msh")
    solid2   = FEM2D("rolling-cyl-plane.msh")

    # as the mesh was created with the center of the disk at (0,0)
	#move(solid1,SVector{2,Float64}([ 0.2+grid.dx  0.2+grid.dx]))
	#move(solid2,SVector{2,Float64}([ 0.8-grid.dx  0.8-grid.dx]))
	Fem.move(solid1,SVector{2,Float64}([ ra+5,  ra+130]))
	Fem.move(solid2,SVector{2,Float64}([ 0.,  2.]))

    Fem.assign_velocity(solid1,[1e3,0.])

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

    fixYForBottom(grid)
    fixYNodes(solid2, "BotSurface")

    solids = [solid1, solid2]
    mats   = [material1,material2]

    Tf       = 0.05
    interval = 1
	dtime    = 0.00005;

	#output1  = PyPlotOutput(interval,"twodisks-results/","Two Disks Collision",(4., 4.))
	output2  = VTKOutput(interval,"rolling-cyl-femp/",["pressure"])
	fix      = EnergiesFix(solids,"rolling-cyl-femp/energies.txt")

    algo2    = USL(1e-9)
    algo1    = TLFEM(1e-9,1.0)
    body     = ConstantBodyForce2D(bodyf)

	report(grid,solids,dtime)

    plotGrid(output2,grid)
    plotParticles_2D(output2,solids,mats,0)

	#reset_timer!
    solve_explicit_dynamics_femp_2D_Contact(grid,solids,mats,basis,body,fric,algo1,output2,fix,Tf,dtime)
    #print_timer()
	# plotting energies


# end

# @time main()