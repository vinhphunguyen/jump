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

module FemMPM

using StaticArrays
using TimerOutputs
using LinearAlgebra
using Printf

using Fix
using Algorithm
using Basis
using Mesh
using Material
using Output

include("SolveExplicitFEMP2D.jl")
include("SolveExplicitFEMP3D.jl")

export solve_explicit_dynamics_femp_2D
export solve_explicit_dynamics_femp_3D

end
