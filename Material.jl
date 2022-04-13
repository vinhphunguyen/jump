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
density::Float64
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

	vmStr    ::Vector{Float64}             # equivalent von Mises stress
	dam      ::Vector{Float64}             # equivalent von Mises stress

	c_dil    :: Float64

	function ElasticMaterial(E,nu,density,Gf,l,parCount)
		lambda = E*nu/(1+nu)/(1-2*nu)
		mu     = E/2/(1+nu)

		c_dil  = sqrt( (lambda+2*mu)/density )

        new(E,nu,density,lambda,mu,Gf,l,0.5/E,fill(0.,parCount),fill(0.,parCount),c_dil)
    end
end

# outer constructor to simplify 
ElasticMaterial(E,nu,rho,parCount)=ElasticMaterial(E,nu,rho,0.,0.,parCount)

# do not update solid.stress!!!
function update_stress!(sigma::MMatrix{2,2,Float64},mat::ElasticMaterial,
	                    epsilon::SMatrix{2,2,Float64},strain_increment,F, J,ip,dtime)
  sigma    .= mat.lambda * (epsilon[1,1]+epsilon[2,2]) * UniformScaling(1.) +
              2.0 * mat.mu * epsilon
end

function update_stress!(sigma::SMatrix{3,3,Float64},mat::ElasticMaterial,
						epsilon::SMatrix{3,3,Float64},strain_increment,F, J,ip,dtime)
  #sigma    .= mat.lambda * (epsilon[1,1]+epsilon[2,2]+epsilon[3,3]) * UniformScaling(1.) + 2.0 * mat.mu * epsilon

  sigma_dev     = sigma - 0.333333333 *(sigma[1,1]+sigma[2,2]+sigma[3,3]) * UniformScaling(1.); 
  mat.vmStr[ip] = sqrt(  sigma[1,1]^2 + sigma[2,2]^2 + sigma[3,3]^2 + 2*(sigma[1,2]^2+sigma[1,3]^2+sigma[2,3]^2) )
  mat.lambda * (epsilon[1,1]+epsilon[2,2]+epsilon[3,3]) * UniformScaling(1.) + 2.0 * mat.mu * epsilon
end

function computeCrackDrivingForce(stress,mat::ElasticMaterial)
	sigma1 = computeMaxPrincipleStress(stress)
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
   return (1.0/J)*( mat.mu*(F*F'-UniformScaling(1.)) + mat.lambda*log(J)*UniformScaling(1.) )
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
	m      :: Float64
	χ      :: Float64  # Taylor-Quinney coeffcient
	Cp     :: Float64  # specific heat
	# Lame's constant
    lambda ::Float64
    mu     ::Float64
    kappa  ::Float64
    signal_velocity::Float64
    cellsize::Float64
    delta::Float64


	alpha    ::Vector{Float64}             # equivalent plastic strain
	alphadot ::Vector{Float64}             # equivalent plastic strain rate
	vmStr    ::Vector{Float64}             # equivalent von Mises stress

	sigmaTrial    ::SMatrix{3,3,Float64,9}
	sigmaTrial_dev::SMatrix{3,3,Float64,9}
	sigmaFinal_dev::SMatrix{3,3,Float64,9}

	c_dil   :: Float64

	function JohnsonCookMaterial(E,nu,density,A,B,C,n,eps0dot,m,χ,Cp,cellsize,parCount)

		lambda = E*nu/(1+nu)/(1-2*nu)
		mu     = E/2/(1+nu)
		kappa  = lambda + mu

		velo   = sqrt((lambda+2*mu)/density)

		alpha    = fill(0,parCount)
		alphad   = fill(0,parCount)
		vmStr    = fill(0,parCount)		

		delta    = χ/(density*Cp)

		c_dil  = sqrt( (lambda+2*mu)/density )

        new(E,nu,density,A,B,C,n,eps0dot,m,χ,Cp,lambda,mu,kappa,velo,cellsize,delta,
        	alpha,alphad,vmStr,zeros(3,3),zeros(3,3),zeros(3,3),c_dil)
    end
end


###############################################################
# VoceMaterial: 
###############################################################
struct VoceMaterial <: MaterialType
    E      ::Float64
    nu     ::Float64
    density::Float64
    A      ::Float64
    B      :: Float64
    C      :: Float64
    χ      :: Float64  # Taylor-Quinney coeffcient
    Cp     :: Float64  # specific heat
    # Lame's constant
    lambda ::Float64
    mu     ::Float64
    kappa  ::Float64
    signal_velocity::Float64
    cellsize::Float64
    delta::Float64


    alpha    ::Vector{Float64}             # equivalent plastic strain
    alphadot ::Vector{Float64}             # equivalent plastic strain rate
    vmStr    ::Vector{Float64}             # equivalent von Mises stress

    sigmaTrial    ::SMatrix{3,3,Float64,9}
    sigmaTrial_dev::SMatrix{3,3,Float64,9}
    sigmaFinal_dev::SMatrix{3,3,Float64,9}

    c_dil   :: Float64

    function VoceMaterial(E,nu,density,A,B,C,χ,Cp,cellsize,parCount)

	lambda = E*nu/(1+nu)/(1-2*nu)
	mu     = E/2/(1+nu)
	kappa  = lambda + mu

	velo   = sqrt((lambda+2*mu)/density)

	alpha    = fill(0,parCount)
	alphad   = fill(0,parCount)
	vmStr    = fill(0,parCount)		

	delta    = χ/(density*Cp)

	c_dil  = sqrt( (lambda+2*mu)/density )

        new(E,nu,density,A,B,C,χ,Cp,lambda,mu,kappa,velo,cellsize,delta,
            alpha,alphad,vmStr,zeros(3,3),zeros(3,3),zeros(3,3),c_dil)
    end
end

# 3D stress update for J2 plasticity
function update_stress!(sig::SMatrix{3,3,Float64},mat::JohnsonCookMaterial,
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

	tav = 1000 * mat.cellsize / mat.signal_velocity;
    mat.alphadot[ip] -= mat.alphadot[ip] * dtime / tav;
    mat.alphadot[ip] += plastic_strain_increment / tav;
    mat.alphadot[ip] = max(0.0, mat.alphadot[ip]);

    mat.vmStr[ip] = sqrt(  sigmaFinal_dev[1,1]^2 + sigmaFinal_dev[2,2]^2 + sigmaFinal_dev[3,3]^2 + 2*(sigmaFinal_dev[1,2]^2+sigmaFinal_dev[1,3]^2+sigmaFinal_dev[2,3]^2) )

    # this is sigma to return
    sigmaFinal_dev     + mat.kappa*(eps0[1,1]+eps0[2,2]+eps0[3,3])*UniformScaling(1.)
end   


# 3D stress update for J2 plasticity
function update_stress!(sig::SMatrix{3,3,Float64},mat::VoceMaterial,
	                eps0::SMatrix{3,3,Float64},strain_increment,F, J, ip,dtime)

    sigmaTrial     = mat.sigmaTrial
    sigmaTrial_dev = mat.sigmaTrial_dev
    sigmaFinal_dev = mat.sigmaFinal_dev

    eff_plastic_strain  = mat.alpha[ip]

    if (eff_plastic_strain < 1.0e-10)
      yieldStress = mat.A;
    else
        yieldStress = mat.A - mat.B * exp(-mat.C * eff_plastic_strain)
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

    tav = 1000 * mat.cellsize / mat.signal_velocity;
    mat.alphadot[ip] -= mat.alphadot[ip] * dtime / tav;
    mat.alphadot[ip] += plastic_strain_increment / tav;
    mat.alphadot[ip] = max(0.0, mat.alphadot[ip]);

    mat.vmStr[ip] = sqrt(  sigmaFinal_dev[1,1]^2 + sigmaFinal_dev[2,2]^2 + sigmaFinal_dev[3,3]^2 + 2*(sigmaFinal_dev[1,2]^2+sigmaFinal_dev[1,3]^2+sigmaFinal_dev[2,3]^2) )

    # this is sigma to return
    sigmaFinal_dev     + mat.kappa*(eps0[1,1]+eps0[2,2]+eps0[3,3])*UniformScaling(1.)
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
	temp    :: Vector{Float64}
	function JohnsonCookDamage(d1,d2,d3,d4,d5,epsdot0,Tr,Tm,parCount)
		return new(d1,d2,d3,d4,d5,epsdot0,Tr,Tm,Tm-Tr,fill(Tr,parCount))
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
    triax = -pH / (vm + 0.01 * abs(pH)); # have softening in denominator to avoid divison by zero
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

  if (mat.d5 > 0.0 && T >= mat.Tr) jc_failure_strain *= 1. + mat.d5 * (T - mat.Tr) / mat.Tmr; end

  damage_init += plastic_strain_increment / jc_failure_strain;

  if (damage_init >= 1.0)
    damage = min((damage_init - 1.0) * 10, 1.0);
  end

  return (damage_init,damage)
end

struct JohnsonCookMaterialWithDamage
	strength_part::JohnsonCookMaterial
	damage_part  ::JohnsonCookDamage


	damage_init    ::Vector{Float64}             
	dam            ::Vector{Float64}     
	density ::Float64
	function   JohnsonCookMaterialWithDamage(strength,damage)
		damage_init = zeros(length(strength.alpha))
		return new(strength,damage,damage_init,copy(damage_init),strength.density)
	end
end

# 3D stress update for J2 plasticity
function update_stress!(sig::MMatrix{3,3,Float64},mat::JohnsonCookMaterialWithDamage,
	                    eps0::SMatrix{3,3,Float64},strain_increment,F, J, ip,dtime)

    damage = mat.dam[ip]

    # id damage = 1: only pressure (no deviatoric stress)
	if (damage >= 1.0)
	    pH  = mat.strength_part.kappa*(eps0[1,1]+eps0[2,2]+eps0[3,3])
	    sig        .=  SMatrix{3,3}(pH,0,0,0,pH,0,0,0,pH)
	    return nothing;
	end


    sigmaTrial     = mat.strength_part.sigmaTrial
    sigmaTrial_dev = mat.strength_part.sigmaTrial_dev
    sigmaFinal_dev = mat.strength_part.sigmaFinal_dev

    eff_plastic_strain  = mat.strength_part.alpha[ip]

    epsdot_ratio        = mat.strength_part.alphadot[ip] / mat.strength_part.epsdot0;
    epsdot_ratio        = max(epsdot_ratio, 1.0);

    T = mat.damage_part.temp[ip]

    if (eff_plastic_strain < 1.0e-10)
      yieldStress = mat.strength_part.A;
    elseif ( T < mat.damage.Tm ) 
      if (mat.strength_part.C != 0)
        yieldStress = (mat.strength_part.A + mat.strength_part.B * eff_plastic_strain^mat.strength_part.n) * (1.0 + epsdot_ratio)^mat.strength_part.C;
      else
        yieldStress = mat.strength_part.A + mat.strength_part.B * eff_plastic_strain^mat.strength_part.n;
      end  
      if (mat.strength_part.m != 0 && T >= mat.damage_part.Tr)
        yieldStress *= 1.0 - ((T - mat.damage_part.Tr) / mat.damage_part.Tmr)^mat.strength_part.m ;
      end
    else
        yieldStress = 0.
    end

    
    Gd     = mat.strength_part.mu
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

	mat.strength_part.alpha[ip]        = eff_plastic_strain + plastic_strain_increment
	# be careful with .=, without it sig is not updated 
	
	pH                   = mat.strength_part.kappa*(eps0[1,1]+eps0[2,2]+eps0[3,3])
	sig                 .= sigmaFinal_dev     + pH*UniformScaling(1.)

	tav = 1000 * mat.strength_part.cellsize / mat.strength_part.signal_velocity;
    mat.strength_part.alphadot[ip] -= mat.strength_part.alphadot[ip] * dtime / tav;
    mat.strength_part.alphadot[ip] += plastic_strain_increment / tav;
    mat.strength_part.alphadot[ip] = max(0.0, mat.strength_part.alphadot[ip]);

    flow_stress            = 1.224744871391589 * sqrt(  sigmaFinal_dev[1,1]^2 + sigmaFinal_dev[2,2]^2 + sigmaFinal_dev[3,3]^2 + 2*(sigmaFinal_dev[1,2]^2+sigmaFinal_dev[1,3]^2+sigmaFinal_dev[2,3]^2) )
    mat.strength_part.vmStr[ip] = flow_stress

    # update damage
    mat.damage_init[ip], mat.dam[ip] = compute_damage(mat.damage_init[ip], mat.dam[ip], pH, sigmaFinal_dev, 
    	                                         mat.strength_part.alphadot[ip], plastic_strain_increment,T, mat.damage_part)

    T += mat.strength_part.delta*flow_stress*plastic_strain_increment;
    mat.damage.temp[ip] = T
end	                    

function getPlasticStrain(ip,mat)
	return mat.alpha[ip]
end

function getPlasticStrain(ip,mat::JohnsonCookMaterialWithDamage)
	return mat.strength.alpha[ip]
end

function getPlasticStrain(ip,mat::Union{ElasticMaterial,RigidMaterial,NeoHookeanMaterial})
	return 0.
end

function getTemperature(ip,mat::Union{ElasticMaterial,RigidMaterial,NeoHookeanMaterial,JohnsonCookMaterial,VoceMaterial})
	return 0.
end

function getTemperature(ip,mat::JohnsonCookMaterialWithDamage)
	return mat.damage.temp[ip]
end


function getDamage(ip,mat::Union{NeoHookeanMaterial,RigidMaterial,JohnsonCookMaterial, VoceMaterial})
	return 0.
end


function getDamage(ip,mat::Union{ElasticMaterial,JohnsonCookMaterialWithDamage})
	return mat.dam[ip]
end

function get_von_mises_stress(ip,mat)
	return mat.vmStr[ip]
end

function get_von_mises_stress(ip,mat::ElasticMaterial)
	return mat.vmStr[ip]
end

function get_von_mises_stress(ip,mat::RigidMaterial)
	return 0.
end


function get_von_mises_stress(ip,mat::JohnsonCookMaterialWithDamage)
	return mat.strength.vmStr[ip]
end



export MaterialType, RigidMaterial, ElasticMaterial, ElastoPlasticMaterial, NeoHookeanMaterial, JohnsonCookMaterial, VoceMaterial, JohnsonCookDamage, JohnsonCookMaterialWithDamage
export update_stress!, computeCrackDrivingForce, getPlasticStrain, get_von_mises_stress,getTemperature, getDamage

end
