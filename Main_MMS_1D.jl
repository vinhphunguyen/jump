# Phu Nguyen, Monash University
# 20 March, 2020 (Coronavirus outbreak)

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/jMPM/src")
# import Gadfly
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")


#include("./Grid.jl")
#include("./Problem.jl")

using Material
using Algorithm
using Solid
using Grid
using Problem
using Output
using Fix
using BodyForce
using Basis
using Mesh

#function main()
	G    = 10e-4

	Rho  = 1000.
    E    = 1e7
    Nu   = 0.3

	c         = sqrt(E/Rho)

	# grid creation and basis
	grid  = Grid1D(1.0, 21)
	basis = LinearBasis()

	fOffset = grid.dx/2

    # do not forget to shift the solid to the right one cell (to have ghost cell)
	coords      = buildParticleForSegment(1., fOffset)
	material    = NeoHookeanMaterial(E,Nu,Rho)
	#println(nodes)
	solid1      = Solid1D(coords,material)

    solid1.mass          .*= fOffset
    solid1.volume        .=  fOffset
    solid1.volumeInitial .=  fOffset

	# assign initial velocity
	for p=1:solid1.parCount
		x = solid1.pos[p]
		solid1.velocity[p] = pi*c*G*sin(pi*x)
	end

    solids = [solid1]

	@printf("	Disk,   number of material points: %d \n", solid1.parCount)

    # Boundary conditions
    fixXForLeft(grid)
    fixXForRight(grid)

    dtime   = 0.1*grid.dx/c;
    Tf      = 2/c
    interval= 20

	bodyforce = AxisAlignBodyForce1D(G,E,Nu,Rho)

	#output1  = PyPlotOutput(interval,"impact-results/","Impact",(4., 4.))
	#output2  = OvitoOutput(interval,"vibratingbeam-cpdi-results/",[])
    problem   = ExplicitSolidMechanics1D(grid,solids,basis,Tf,bodyforce,[])
    algo      =  MUSL()
    solve(problem, algo, interval, dtime)
    ############################################################
    # plot solutions
	############################################################
    step100 = 3
	t= problem.recordTime[step100]
	u_exact = Vector{Float64}(undef,0)
	xx = solid1.pos0
    for p=1:solid1.parCount
        uu = G*sin(pi*xx[p])*sin(c*pi*t)
		push!(u_exact,uu)
	end

	pyFig_RealTime = PyPlot.figure("MPM 2Disk FinalPlot", figsize=(8/2.54, 4/2.54))
	PyPlot.clf()
	pyPlot01 = PyPlot.gca()
	PyPlot.subplots_adjust(left=0.15, bottom=0.25, right=0.65)
	pyPlot01[:grid](b=true, which="both", color="gray", linestyle="-", linewidth=0.5)
	pyPlot01[:set_axisbelow](true)
	# pyPlot01[:set_xlim](0.0, 1.0)
	# pyPlot01[:set_ylim](0.0, 1.2*G)
	# pyPlot01[:set_xlabel]("position", fontsize=8)
	# pyPlot01[:set_ylabel]("displacement", fontsize=8)
	# pyPlot01[:set_xticks](collect(0.0:1.0))
	# pyPlot01[:tick_params](axis="both", which="major", labelsize=8)
	# pyPlot01[:set_yticks](collect(0.0:G))
	PyPlot.plot(xx, c="blue", u_exact, "-", label="\$ K \$", linewidth=1.0)
	PyPlot.plot(xx, c="red", problem.recordX[:,step100]-solid1.pos0, ".", label="\$ K \$", linewidth=1.0)
	#PyPlot.plot(xx, c="cyan", problem.recordX[:,3]-solid1.pos0, "-", label="\$ K \$", linewidth=1.0)
	PyPlot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=8)
	#
	PyPlot.savefig("plot_2Disk_Julia.pdf")
# end
#
# @time main()
