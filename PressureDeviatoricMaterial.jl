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

abstract type EOSType       end
abstract type StrengthType  end

struct PressureDeviatoricMaterial{P<:EOSType,D:StrengthType}
	pressure    :: P
	dev_stress  :: D
	function PressureDeviatoricMaterial(p::P,d::D) where {P<:EOSType,D:StrengthType}
		return PressureDeviatoricMaterial{P,D}(pressure,dev_stress)
	end
end

function update_stress!(ip, solid, vel_grad, dtime, mat::PressureDeviatoricMaterial{P,D}) where {P<:EOSType,D:StrengthType}
      #mat->eos->compute_pressure(pH, ienergy[ip], J[ip], rho[ip], T[ip],damage[ip]);
      # sigma_dev = mat->strength->update_deviatoric_stress(
      #     sigma[ip], D[ip], plastic_strain_increment, eff_plastic_strain[ip],
      #     eff_plastic_strain_rate[ip], damage[ip], T[ip]);

      # eff_plastic_strain[ip] += plastic_strain_increment;

      # // // compute a characteristic time over which to average the plastic
      # // strain
      # double tav = 1000 * grid->cellsize / mat->signal_velocity;
      # eff_plastic_strain_rate[ip] -=
      #     eff_plastic_strain_rate[ip] * update->dt / tav;
      # eff_plastic_strain_rate[ip] += plastic_strain_increment / tav;
      # eff_plastic_strain_rate[ip] = MAX(0.0, eff_plastic_strain_rate[ip]);

      # if (mat->damage != NULL)
      #   mat->damage->compute_damage(damage_init[ip], damage[ip], pH, sigma_dev,
      #                               eff_plastic_strain_rate[ip],
      #                               plastic_strain_increment, T[ip]);
      # sigma[ip] = -pH * eye + sigma_dev;	

      update_pressure!(p)
      update_deviatoric_stress!(sigma_dev)

      solid.stress[ip] = -p*UniformScaling(1.) + sigma_dev;

end

struct LinearEOS <: EOSType
	rho::Float64
	K  ::FLoat64
end

struct FluidEOS <: EOSType
	rho::Float64
	K  ::FLoat64
	γ  ::FLoat64
end

struct FluidStrength <: StrengthType
	G :: Float64
end

struct JohnsonCook <: StrengthType
  G       :: Float64
  A       :: Float64
  B       :: Float64
  n       :: Float64
  epsdot0 :: Float64
  C       :: Float64
  m       :: Float64
  Tr      :: Float64
  Tm      :: Float64
end

function update_pressure!(,:mat:LinearEOS)
    return  K*(1-J)
end

function update_pressure!(,mat::FluidEOS)
   mu = rho / mat.rho
   return mat.K * (pow(mu, mat.γ) - 1.0)
end

function update_deviatoric_stress!(sigma_dev,mat::FluidStrength)

end
#=
function update_deviatoric_stress!(sigma_dev,mat::JohnsonCook)
   A = mat.A
   B = mat.B
   C = mat.C

Matrix3d StrengthJohnsonCook::update_deviatoric_stress(
    const Eigen::Matrix3d &sigma, const Eigen::Matrix3d &D,
    double &plastic_strain_increment, const double eff_plastic_strain,
    const double epsdot, const double damage, const double T)
  Matrix3d sigmaFinal_dev, sigmaTrial, sigmaTrial_dev;
  double J2, Gd, yieldStress;

  double epsdot_ratio = epsdot / epsdot0;
  epsdot_ratio        = MAX(epsdot_ratio, 1.0);

  if (eff_plastic_strain < 1.0e-10) 
    yieldStress = A
  else 
    if (T < mat.Tm)
      if (C != 0)
        yieldStress = (A + B * pow(eff_plastic_strain, mat.n)) * pow(1.0 + epsdot_ratio, C);
      else
        yieldStress = A + B * pow(eff_plastic_strain, mat.n);
      end
      if (m != 0 && T >= Tr)
        yieldStress *= 1.0 - pow((T - Tr) / Tmr, m);
    else
      yieldStress = 0
  end
  end

  /*
   * deviatoric rate of unrotated stress
   */
  Gd               = G_;

  // sigmaInitial_dev = Deviator(sigma);

  /*
   * perform a trial elastic update to the deviatoric stress
   */
  sigmaTrial = sigma + update->dt * 2.0 * Gd * D;
  sigmaTrial_dev = Deviator(sigmaTrial); // increment stress deviator using deviatoric rate

  /*
   * check yield condition
   */
  J2             = SQRT_3_OVER_2 * sigmaTrial_dev.norm();
  sigmaFinal_dev = sigmaTrial_dev;

  if (J2 < yieldStress)
  {
    /*
     * no yielding has occured.
     * final deviatoric stress is trial deviatoric stress
     */
    plastic_strain_increment = 0.0;
  }
  else
  {
    // printf("yiedl\n");
    /*
     * yielding has occured
     */

    plastic_strain_increment = (J2 - yieldStress) / (3.0 * Gd);
    /*
     * new deviatoric stress:
     * obtain by scaling the trial stress deviator
     */
    sigmaFinal_dev *= (yieldStress / J2);
  }
=#  

export FluidEOS, LinearEOS
export update_stress!

end
