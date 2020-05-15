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
module PressureDeviatoricMaterial

using TimerOutputs
using StaticArrays
using LinearAlgebra


using Util

abstract type PressureType         end
abstract type DeviatoricStressType end

struct PressureDeviatoricMaterial{P<:PressureType,D:DeviatoricStressType}
	pressure    :: P
	dev_stress  :: D
	function PressureDeviatoricMaterial(p::P,d::D) where {P<:PressureType,D:DeviatoricStressType}
		return PressureDeviatoricMaterial{P,D}(pressure,dev_stress)
	end
end

function update_stress!(ip, solid, vel_grad, dtime, mat)
	D                  = 0.5 * (vel_grad + vel_grad')
	solid.strain[ip]  += dtime * D
	F[ip]             *= (UniformScaling(1.) + vel_grad*dtime)
	J                  = det(F[ip])
	solid.volume[ip]   = J * solid.volumeInitial[ip]
	epsilon            = solid.strain[ip]
 	
    soldi.stress[ip]    .= mat.lambda * (epsilon[1,1] + epsilon[2,2]) * UniformScaling(1.) +
                        2.0 * mat.mu * epsilon
end


end
