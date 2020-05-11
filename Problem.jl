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

# Module Problem
#

module Problem

using PyPlot
using StaticArrays
using LinearAlgebra
using Printf

using Grid
using Solid
using Basis
using Material
using Output
using Algorithm
using Fix
using BodyForce

# #########################################################
# #  1D problem
# #########################################################
#
# mutable struct ExplicitSolidMechanics1D
# 	grid      ::Grid1D
# 	solids    ::Vector{Solid1D}
# 	basis     :: BasisType
# 	solidCount::Int64
# 	Tf        ::Float64
# 	bodyforce ::BodyForceType
# 	fixes     ::FixBase
# 	kinEnergy ::Vector{Float64}
# 	strEnergy ::Vector{Float64}
# 	recordTime::Vector{Float64}
# 	recordX   ::Array{Float64,2}
#
#   function ExplicitSolidMechanics1D(grid::Grid1D,solids::Vector{Solid1D}, basis,
# 	                              Tf::Float64,bodyforce,fixes)
#   	new(grid,solids,basis,length(solids),Tf,bodyforce,fixes,
# 	   Vector{Float64}(),Vector{Float64}(),Vector{Float64}(),Array{Float64,2}(undef,solids[1].parCount,0))
#   end
# end
#
#
# #########################################################
# #  2D problem, ULMPM
# #########################################################
#
# struct ExplicitSolidMechanics2D{T <: BasisType}
# 	grid      ::Grid2D
# 	solids    ::Vector{Solid2D}
# 	basis     :: T
# 	solidCount::Int64
# 	Tf        ::Float64
# 	#bodyforce ::BodyForceType
# 	output    ::OutputType
# 	fixes     ::FixBase      # Vector{EmptyFix} cannot be used for Vector{FixBase}
# 	kinEnergy ::Vector{Float64}
# 	strEnergy ::Vector{Float64}
# 	recordTime::Vector{Float64}
#
#   function ExplicitSolidMechanics2D(grid::Grid2D,solids::Vector{Solid2D}, basis::T,
# 	               Tf::Float64,output::OutputType,fixes::FixBase) where {T<:BasisType}
#   	return new{T}(grid,solids,basis,length(solids),Tf,output,fixes,
# 	    Vector{Float64}(),Vector{Float64}(),Vector{Float64}())
#   end
# end
#
# #########################################################
# #  2D problem, TLMPM
# #########################################################
#
# struct ExplicitDynamicsTL2D
# 	grids     ::Vector{Grid2D}
# 	solids    ::Vector{Solid2D}
# 	basis     :: BasisType
# 	solidCount::Int64
# 	Tf        ::Float64
# 	bodyforce ::BodyForceType
# 	output    ::OutputType
# 	fixes     ::FixBase
# 	kinEnergy ::Vector{Float64}
# 	strEnergy ::Vector{Float64}
# 	recordTime::Vector{Float64}
#
#   function ExplicitDynamicsTL2D(grids::Vector{Grid2D},solids::Vector{Solid2D}, basis,
# 	                              Tf::Float64,bodyforce,output::OutputType,fixes)
#   	new(grids,solids,basis,length(solids),Tf,bodyforce,output,fixes,Vector{Float64}(),Vector{Float64}(),Vector{Float64}())
#   end
# end
#
# #########################################################
# #  3D problem
# #########################################################
#
# struct ExplicitSolidMechanics3D
# 	grid      ::Grid3D
# 	solids    ::Vector{Solid3D}
# 	basis     :: BasisType
# 	solidCount::Int64
# 	Tf        ::Float64
# 	bodyforce ::BodyForceType
# 	output    ::OutputType
#     fixes
#
# 	kinEnergy ::Vector{Float64}
# 	strEnergy ::Vector{Float64}
# 	recordTime::Vector{Float64}
#
#   function ExplicitSolidMechanics3D(grid::Grid3D,solids::Vector{Solid3D},basis, Tf::Float64,
# 	                                bodyforce, output::OutputType, fixes )
#   	new(grid,solids,basis,length(solids),Tf,bodyforce,output,fixes,Vector{Float64}(),Vector{Float64}(),Vector{Float64}())
#   end
# end
#
# #########################################################
# #  2D problem, TLMPM
# #########################################################
#
# struct ExplicitDynamicsPhaseField2D
# 	grid      ::Grid2D
# 	solids    ::Vector{Solid2D}
# 	basis     :: BasisType
# 	solidCount::Int64
# 	Tf        ::Float64
# 	bodyforce ::BodyForceType
# 	output    ::OutputType
# 	fixes     ::FixBase
# 	kinEnergy ::Vector{Float64}
# 	strEnergy ::Vector{Float64}
# 	recordTime::Vector{Float64}
#
#   function ExplicitDynamicsPhaseField2D(grid::Grid2D,solids::Vector{Solid2D}, basis,
# 	                              Tf::Float64,bodyforce,output::OutputType,fixes)
#   	new(grid,solids,basis,length(solids),Tf,bodyforce,output,fixes,Vector{Float64}(),Vector{Float64}(),Vector{Float64}())
#   end
# end

# these files implement the solve(Problem,....)
# for 1D, 2D and 3D problems
include("SolveExplicitSolidMechanics1D.jl")
include("SolveExplicitSolidMechanics2D.jl")
include("SolveExplicitSolidMechanics3D.jl")

include("SolveExplicitDynamicsTL2D.jl")
#include("SolveExplicitDynamicsTL3D.jl")

include("SolveExplicitDynamicsPFM2D.jl")

export solve_explicit_dynamics_2D

end
