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
using Fem

include("SolveExplicitFEMP2D.jl")
include("SolveExplicitFEMP3D.jl")
include("SolveExplicitFEMP2DContact.jl")

export solve_explicit_dynamics_femp_2D
export solve_explicit_dynamics_femp_3D

export solve_explicit_dynamics_femp_2D_Contact
#export solve_explicit_dynamics_femp_3D_Contact

end
