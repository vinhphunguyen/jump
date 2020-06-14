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
module Material

using TimerOutputs
using StaticArrays
using LinearAlgebra


using Util

abstract type MaterialType end

Identity   = UniformScaling(1.)


###############################################################
# RigidMaterial: model rigid bodies
###############################################################
struct RigidMaterial <: MaterialType
	vx
	vy
	density    # not used, but needed for compability with density of non-rigid materials
	fixed::Bool
	function RigidMaterial(vx,vy,density)
		fixed = false
		if vx == 0. && vy == 0. fixed = true end
		return new(vx,vy,density,fixed)
	end
end


###############################################################
# ElasticMaterial: Hooke materials
###############################################################
struct ElasticMaterial <: MaterialType
	E      ::Float64
	nu     ::Float64
	density::Float64
    # Lame's constant
    lambda ::Float64
    mu     ::Float64
	Gf     ::Float64
	l0     ::Float64
	fac    ::Float64

	function ElasticMaterial(E,nu,density,Gf,l)
		lambda = E*nu/(1+nu)/(1-2*nu)
		mu     = E/2/(1+nu)

        new(E,nu,density,lambda,mu,Gf,l,0.5/E)
    end
end

# outer constructor to simplify 
ElasticMaterial(E,nu,rho)=ElasticMaterial(E,nu,rho,0.,0.)

# do not update solid.stress!!!
function update_stress!(sigma::MMatrix{2,2,Float64},mat::ElasticMaterial,
	                    epsilon::SMatrix{2,2,Float64},strain_increment,F, J,ip,dtime)
  sigma    .= mat.lambda * (epsilon[1,1]+epsilon[2,2]) * UniformScaling(1.) +
              2.0 * mat.mu * epsilon
end

function update_stress!(sigma::MMatrix{3,3,Float64},mat::ElasticMaterial,
						epsilon::SMatrix{3,3,Float64},strain_increment,F, J,ip,dtime)
  sigma    .= mat.lambda * (epsilon[1,1]+epsilon[2,2]+epsilon[3,3]) * UniformScaling(1.) + 2.0 * mat.mu * epsilon
end

function computeCrackDrivingForce(stress,mat::ElasticMaterial)
	sigma1 = computeMaxPrincipleStress2D(stress)
	term1  = max( sigma1, 0. )

    return mat.fac * term1 * term1;
end

# function update_stress(mat::ElasticMaterial,epsilon::SMatrix{2,2,Float64})
#      return mat.lambda * (epsilon[1,1]+epsilon[2,2]) * Identity + 2.0 * mat.mu * epsilon
#    end


###############################################################
# NeoHookeanMaterial: compressible Neo-Hookean Material
###############################################################
 struct NeoHookeanMaterial <: MaterialType
 	E      ::Float64
 	nu     ::Float64
 	density::Float64
 	# Lame's constant
 	lambda ::Float64
 	mu     ::Float64

 	function NeoHookeanMaterial(E,nu,density)

 		lambda = E*nu/(1+nu)/(1-2*nu)
 		mu     = E/2/(1+nu)

 		new(E,nu,density,lambda,mu)
 	end
 end
 # do not update solid.stress!!!
 function update_stress!(sigma,mat::NeoHookeanMaterial,
	                     epsilon, strain_increment, F, J, ip,dtime)
   sigma    .= (1.0/J)*( mat.mu*(F*F'-UniformScaling(1.)) + mat.lambda*log(J)*UniformScaling(1.) )
 end


###############################################################
# ElastoPlasticMaterial: 2D J2 plastic materials (metals,...)
###############################################################
struct ElastoPlasticMaterial <: MaterialType
  	E      ::Float64
	nu     ::Float64
	density::Float64
	# Lame's constant
    lambda ::Float64
    mu     ::Float64
    kappa  ::Float64
    # plastic constants
	fy     ::Float64         # yield stress
	k1     ::Float64         # hardening modulus  SMatrix{2,2,Float64,4}

	pstrain::Vector{SVector{3,Float64}}  # plastic strain
	alpha  ::Vector{Float64}             # equivalent plastic strain
	vmStr  ::Vector{Float64}             # equivalent von Mises stress

    # 2D version
	epsilon_dev::MVector{3,Float64}
	s_trial    ::MVector{3,Float64}
	sigma_trial::MVector{3,Float64}
	dev_Sig    ::MVector{3,Float64}

    # 3D version
    pstrain_3D::Vector{MMatrix{3,3,Float64,9}}  # plastic strain
	epsilon_dev_3D::MMatrix{3,3,Float64,9}
	s_trial_3D    ::MMatrix{3,3,Float64,9}
	sigma_trial_3D::MMatrix{3,3,Float64,9}
	dev_Sig_3D    ::MMatrix{3,3,Float64,9}

	function ElastoPlasticMaterial(E,nu,density,fy,k1,parCount)

		lambda = E*nu/(1+nu)/(1-2*nu)
		mu     = E/2/(1+nu)
		kappa  = lambda + mu

		alpha    = fill(0,parCount)
		vmStr    = fill(0,parCount)
		pstrain  = fill(zeros(3),parCount)
		pstrain_3D  = fill(zeros(3,3),parCount)

        new(E,nu,density,lambda,mu,kappa,fy,k1,pstrain,alpha,vmStr,
		zeros(3),zeros(3),zeros(3),zeros(3),pstrain_3D,zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3))
    end
end

# 2D stress update for J2 plasticity
function update_stress!(sig::MMatrix{2,2,Float64},mat::ElastoPlasticMaterial,
	                    eps0::SMatrix{2,2,Float64},strain_increment,F, J, ip,dtime)

	eye2   = @SVector [1., 1., 0.];
	eye2x2 = SMatrix{3,3}(1., 1., 0., 1., 1., 0., 0., 0., 0.)
	I_dev  = SMatrix{3,3}(1., 0., 0., 0., 1., 0., 0., 0., 1.) - 0.5*eye2x2;
	I      = SMatrix{3,3}(1., 0., 0., 0., 1., 0., 0., 0., 0.5) #0.5 to make engineering strain to physical one
	Iinv   = SMatrix{3,3}(1., 0., 0., 0., 1., 0., 0., 0., 2.)

	# Compute trial stress
	epsilonp0    = mat.pstrain[ip]
	alpha0       = mat.alpha[ip]
	epsilon_dev  = mat.epsilon_dev
	s_trial      = mat.s_trial
	sigma_trial  = mat.sigma_trial
	devSig       = mat.dev_Sig

	epsilon      = @SVector [eps0[1,1], eps0[2,2], 2*eps0[1,2]]
	mul!(epsilon_dev,I_dev,epsilon)
	#s_trial      = 2*mat.mu*I*(epsilon_dev-epsilonp0)
	mul!(s_trial,I,@SVector[2*mat.mu*(epsilon_dev[1]-epsilonp0[1]),
	                        2*mat.mu*(epsilon_dev[2]-epsilonp0[2]),
							2*mat.mu*(epsilon_dev[3]-epsilonp0[3])])
	norm_s_trial = sqrt(s_trial[1]^2 + s_trial[2]^2 + 2*s_trial[3]^2)
	sigma_trial  = mat.kappa*(epsilon[1]+epsilon[2])*eye2 + s_trial

	# Check yield condition
	f_trial = norm_s_trial - (mat.k1*alpha0 + mat.fy);

	if f_trial <= 0. # elastic step
	    # alpha    = alpha0;
	    # epsilonp = epsilonp0;
	    sigma    = sigma_trial;
	else # plastic step
	    normal = (1.0/norm_s_trial)* @SVector[s_trial[1], s_trial[2], s_trial[3]]#s_trial/norm_s_trial;
	    lambda = (norm_s_trial - mat.k1*alpha0 - mat.fy)/(2*mat.mu + mat.k1)
	    alpha0 += lambda
	    # Update plastic strain and stress
	    epsilonp0 +=  lambda*@SVector[Iinv[1,1]*normal[1]+Iinv[1,2]*normal[2]+Iinv[1,3]*normal[3],
		                              Iinv[2,1]*normal[1]+Iinv[2,2]*normal[2]+Iinv[2,3]*normal[3],
		                              Iinv[3,1]*normal[1]+Iinv[3,2]*normal[2]+Iinv[3,3]*normal[3]]
	    sigma  = mat.kappa*(epsilon[1]+epsilon[2])*eye2 + s_trial - 2*mat.mu*lambda*normal
	end
	sig[1,1] = sigma[1]
	sig[2,2] = sigma[2]
	sig[1,2] = sigma[3]
	sig[2,1] = sigma[3]
	mat.alpha[ip]    = alpha0
	mat.pstrain[ip]  = epsilonp0
	mul!(devSig,I_dev,sigma)
    mat.vmStr[ip] = sqrt(  devSig[1]^2 + devSig[2]^2 + 2*devSig[3]^2  )
end

# 3D stress update for J2 plasticity
function update_stress!(sig::MMatrix{3,3,Float64},mat::ElastoPlasticMaterial,
	                    eps0::SMatrix{3,3,Float64},strain_increment,F, J, ip,dtime)

	eye2   = @SVector [1., 1., 0.];
	eye2x2 = SMatrix{3,3}(1., 1., 0., 1., 1., 0., 0., 0., 0.)
	I_dev  = SMatrix{3,3}(1., 0., 0., 0., 1., 0., 0., 0., 1.) - 0.5*eye2x2;
	I      = SMatrix{3,3}(1., 0., 0., 0., 1., 0., 0., 0., 0.5) #0.5 to make engineering strain to physical one
	Iinv   = SMatrix{3,3}(1., 0., 0., 0., 1., 0., 0., 0., 2.)

	# Compute trial stress
	epsilonp0    = mat.pstrain_3D[ip]
	alpha0       = mat.alpha[ip]
	epsilon_dev  = mat.epsilon_dev_3D
	s_trial      = mat.s_trial_3D
	sigma_trial  = mat.sigma_trial_3D
	devSig       = mat.dev_Sig_3D

    # deviatoric strain/stress	
	epsilon_dev  = eps0 - 0.333333333 * (eps0[1,1]+eps0[2,2]+eps0[3,3]) * UniformScaling(1.)
	s_trial      = 2*mat.mu*(epsilon_dev-epsilonp0)
	norm_s_trial = sqrt(s_trial[1,1]^2 + s_trial[2,2]^2 + s_trial[3,3]^2+2.0 * (s_trial[1,2]^2 + s_trial[2,3]^2+s_trial[1,3]^2) )
	sigma_trial  = mat.kappa*(eps0[1,1]+eps0[2,2]+eps0[3,3])*UniformScaling(1.) + s_trial


	# Check yield condition
	f_trial = norm_s_trial - (mat.k1*alpha0 + mat.fy);

	if f_trial <= 0. # elastic step
	    # alpha    = alpha0;
	    # epsilonp = epsilonp0;
	    sigma    = sigma_trial;
	else # plastic step
	    normal = (1.0/norm_s_trial) * s_trial
	    lambda = (norm_s_trial - mat.k1*alpha0 - mat.fy)/(2*mat.mu + mat.k1)
	    alpha0 += lambda
	    # Update plastic strain and stress
	    epsilonp0 +=  lambda*normal
	    sigma      = sigma_trial - 2*mat.mu*lambda*normal
	end

	mat.alpha[ip]       = alpha0
	mat.pstrain_3D[ip]  = epsilonp0

	sig .= sigma
	
	devSig   = sigma - 0.333333333 * (sigma[1,1]+sigma[2,2]+sigma[3,3]) * UniformScaling(1.)
    mat.vmStr[ip] = sqrt(  devSig[1,1]^2 + devSig[2,2]^2 + devSig[3,3]^2 + 2*(devSig[1,2]^2+devSig[1,3]^2+devSig[2,3]^2) )
end

###############################################################
# JohnsonCookMaterial: 
###############################################################
struct JohnsonCookMaterial <: MaterialType
  	E      ::Float64
	nu     ::Float64
	density::Float64
	A      ::Float64
	B      :: Float64
	C      :: Float64
	n      :: Float64
	epsdot0:: Float64
	# Lame's constant
    lambda ::Float64
    mu     ::Float64
    kappa  ::Float64
    signal_velocity::Float64
    cellsize::Float64


	alpha    ::Vector{Float64}             # equivalent plastic strain
	alphadot ::Vector{Float64}             # equivalent plastic strain rate
	vmStr    ::Vector{Float64}             # equivalent von Mises stress

	sigmaTrial    ::MMatrix{3,3,Float64,9}
	sigmaTrial_dev::MMatrix{3,3,Float64,9}
	sigmaFinal_dev::MMatrix{3,3,Float64,9}

	function JohnsonCookMaterial(E,nu,density,A,B,C,n,eps0dot,cellsize,parCount)

		lambda = E*nu/(1+nu)/(1-2*nu)
		mu     = E/2/(1+nu)
		kappa  = lambda + mu

		velo   = sqrt((lambda+2*mu)/density)

		alpha    = fill(0,parCount)
		alphad   = fill(0,parCount)
		vmStr    = fill(0,parCount)		

        new(E,nu,density,A,B,C,n,eps0dot,lambda,mu,kappa,velo,cellsize,alpha,alphad,vmStr,zeros(3,3),zeros(3,3),zeros(3,3))
    end
end


# 3D stress update for J2 plasticity
function update_stress!(sig::MMatrix{3,3,Float64},mat::JohnsonCookMaterial,
	                    eps0::SMatrix{3,3,Float64},strain_increment,F, J, ip,dtime)

    sigmaTrial     = mat.sigmaTrial
    sigmaTrial_dev = mat.sigmaTrial_dev
    sigmaFinal_dev = mat.sigmaFinal_dev

    eff_plastic_strain  = mat.alpha[ip]

    epsdot_ratio        = mat.alphadot[ip] / mat.epsdot0;
    epsdot_ratio        = max(epsdot_ratio, 1.0);

    if (eff_plastic_strain < 1.0e-10)
      yieldStress = mat.A;
    else 
      if (mat.C != 0)
        yieldStress = (mat.A + mat.B * eff_plastic_strain^mat.n) * (1.0 + epsdot_ratio)^mat.C;
      else
        yieldStress = mat.A + mat.B * eff_plastic_strain^mat.n;
      end  
    end

    #damage = mat.damage[ip]
    Gd     = mat.mu
	# if (damage > 0)
	# 	Gd          *= (1 - damage);
	# 	yieldStress *= (1 - damage);
	# end

    sigmaTrial     = sig + 2.0 * Gd * strain_increment;
    sigmaTrial_dev = sigmaTrial - 0.333333333 *(sigmaTrial[1,1]+sigmaTrial[2,2]+sigmaTrial[3,3]) * UniformScaling(1.); 

    J2             = 1.224744871391589 * sqrt( sigmaTrial_dev[1,1]^2 + sigmaTrial_dev[2,2]^2 + sigmaTrial_dev[3,3]^2 +
    	                                 2* (sigmaTrial_dev[1,2]^2+sigmaTrial_dev[1,3]^2+sigmaTrial_dev[2,3]^2) )
    sigmaFinal_dev = sigmaTrial_dev;

    if (J2 < yieldStress)
       plastic_strain_increment = 0.0;
    else
       plastic_strain_increment = (J2 - yieldStress) / (3.0 * Gd);
       sigmaFinal_dev     *= (yieldStress / J2);
    end

	mat.alpha[ip]        = eff_plastic_strain + plastic_strain_increment
	# be careful with .=, without it sig is not updated 
	sig                 .= sigmaFinal_dev     + mat.kappa*(eps0[1,1]+eps0[2,2]+eps0[3,3])*UniformScaling(1.)

	tav = 1000 * mat.cellsize / mat.signal_velocity;
    mat.alphadot[ip] -= mat.alphadot[ip] * dtime / tav;
    mat.alphadot[ip] += plastic_strain_increment / tav;
    mat.alphadot[ip] = max(0.0, mat.alphadot[ip]);

    mat.vmStr[ip] = sqrt(  sigmaFinal_dev[1,1]^2 + sigmaFinal_dev[2,2]^2 + sigmaFinal_dev[3,3]^2 + 2*(sigmaFinal_dev[1,2]^2+sigmaFinal_dev[1,3]^2+sigmaFinal_dev[2,3]^2) )
end   

struct JohnsonCookDamage
	d1      :: Float64
	d2      :: Float64
	d3      :: Float64
	d4      :: Float64
	d5      :: Float64
	epsdot0 :: Float64
	Tr      :: Float64
	Tm      :: Float64
	Tmr     :: Float64
	function JohnsonCookDamage(d1,d2,d3,d4,d5,epsdot0,Tr,Tm)
		return new(d1,d2,d3,d4,d5,epsdot0,Tr,Tm,Tm-Tr)
	end
end

function compute_damage(damage_init, damage,
                        pH, Sdev, epsdot, plastic_strain_increment,T, mat::JohnsonCookDamage)

  vm = 1.224744871 *  sqrt( Sdev[1,1]^2 + Sdev[2,2]^2 + Sdev[3,3]^2 + 2* (Sdev[1,2]^2+Sdev[1,3]^2+Sdev[2,3]^2) )

  if (vm < 0.0)
    error("negative von mises stress!!!");
  end

  # determine stress triaxiality
  triax = 0.0;
  if (pH != 0.0 && vm != 0.0)
    triax = -pH / (vm + 0.01 * fabs(pH)); # have softening in denominator to avoid divison by zero
  end

  if (triax > 3.0) triax = 3.0; end
  
  # Johnson-Cook failure strain, dependence on stress triaxiality
  jc_failure_strain = mat.d1 + mat.d2 * exp(mat.d3 * triax);

  # include strain rate dependency if parameter d4 is defined and current
  # plastic strain rate exceeds reference strain rate
  if (mat.d4 > 0.0)
    if (epsdot > mat.epsdot0)
      epdot_ratio        = epsdot / mat.epsdot0;
      jc_failure_strain *= (1.0 + mat.d4 * log(epdot_ratio));
    end
  end

  if (mat.d5 > 0.0 && T >= mat.Tr) jc_failure_strain *= 1 + mat.d5 * (T - mat.Tr) / mat.Tmr; end

  damage_init += plastic_strain_increment / jc_failure_strain;

  if (damage_init >= 1.0)
    damage = min((damage_init - 1.0) * 10, 1.0);
  end

  return (damage_init,damage)
end

struct JohnsonCookMaterialWithDamage
	strength::JohnsonCookMaterial
	damage  ::JohnsonCookDamage


	damage_init    ::Vector{Float64}             
	dam            ::Vector{Float64}     
	function   JohnsonCookMaterialWithDamage(strength,damage)
		damage_init = zeros(length(strength.alpha))
		return new(strength,damage,damage_init,copy(damage_init))
	end
end

# 3D stress update for J2 plasticity
function update_stress!(sig::MMatrix{3,3,Float64},mat::JohnsonCookMaterialWithDamage,
	                    eps0::SMatrix{3,3,Float64},strain_increment,F, J, ip,dtime)
    sigmaTrial     = mat.strength.sigmaTrial
    sigmaTrial_dev = mat.strength.sigmaTrial_dev
    sigmaFinal_dev = mat.strength.sigmaFinal_dev

    eff_plastic_strain  = mat.strength.alpha[ip]

    epsdot_ratio        = mat.strength.alphadot[ip] / mat.strength.epsdot0;
    epsdot_ratio        = max(epsdot_ratio, 1.0);

    if (eff_plastic_strain < 1.0e-10)
      yieldStress = mat.strength.A;
    else 
      if (mat.strength.C != 0)
        yieldStress = (mat.strength.A + mat.strength.B * eff_plastic_strain^mat.strength.n) * (1.0 + epsdot_ratio)^mat.strength.C;
      else
        yieldStress = mat.strength.A + mat.strength.B * eff_plastic_strain^mat.strength.n;
      end  
    end

    damage = mat.damage.dam[ip]
    Gd     = mat.strength.mu
	if (damage > 0)
		Gd          *= (1 - damage);
		yieldStress *= (1 - damage);
	end

    sigmaTrial     = sig + 2.0 * Gd * strain_increment;
    sigmaTrial_dev = sigmaTrial - 0.333333333 *(sigmaTrial[1,1]+sigmaTrial[2,2]+sigmaTrial[3,3]) * UniformScaling(1.); 

    J2             = 1.224744871391589 * sqrt( sigmaTrial_dev[1,1]^2 + sigmaTrial_dev[2,2]^2 + sigmaTrial_dev[3,3]^2 +
    	                                 2* (sigmaTrial_dev[1,2]^2+sigmaTrial_dev[1,3]^2+sigmaTrial_dev[2,3]^2) )
    sigmaFinal_dev = sigmaTrial_dev;

    if (J2 < yieldStress)
       plastic_strain_increment = 0.0;
    else
       plastic_strain_increment = (J2 - yieldStress) / (3.0 * Gd);
       sigmaFinal_dev     *= (yieldStress / J2);
    end

	mat.alpha[ip]        = eff_plastic_strain + plastic_strain_increment
	# be careful with .=, without it sig is not updated 
	
	pH                   = mat.kappa*(eps0[1,1]+eps0[2,2]+eps0[3,3])*UniformScaling(1.)
	sig                 .= sigmaFinal_dev     + pH

	tav = 1000 * mat.cellsize / mat.signal_velocity;
    mat.strength.alphadot[ip] -= mat.alphadot[ip] * dtime / tav;
    mat.strength.alphadot[ip] += plastic_strain_increment / tav;
    mat.strength.alphadot[ip] = max(0.0, mat.alphadot[ip]);

    mat.strength.vmStr[ip] = sqrt(  sigmaFinal_dev[1,1]^2 + sigmaFinal_dev[2,2]^2 + sigmaFinal_dev[3,3]^2 + 2*(sigmaFinal_dev[1,2]^2+sigmaFinal_dev[1,3]^2+sigmaFinal_dev[2,3]^2) )

    # update damage
    damage_init[ip], dam[ip] = compute_damage(damage_init[ip], dam[ip], pH, sigmaFinal_dev, 
    	                                         mat.strength.alphadot[ip], plastic_strain_increment,T, mat.damage)
end	                    

export MaterialType, RigidMaterial, ElasticMaterial, ElastoPlasticMaterial, NeoHookeanMaterial, JohnsonCookMaterial
export update_stress!, computeCrackDrivingForce

end
