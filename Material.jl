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
	vel::Vector{Float64}
	density    # not used, but needed for compability with density of non-rigid materials
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
# do not update solid.stress!!!
function update_stress!(sigma::MMatrix{2,2,Float64},mat::ElasticMaterial,
	                    epsilon::SMatrix{2,2,Float64},F, J,ip)
  sigma    .= mat.lambda * (epsilon[1,1]+epsilon[2,2]) * UniformScaling(1.) +
              2.0 * mat.mu * epsilon
end

function update_stress!(sigma::MMatrix{3,3,Float64},mat::ElasticMaterial,
						epsilon::SMatrix{3,3,Float64},F, J,ip)
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
 function update_stress!(sigma::MMatrix{2,2,Float64},mat::NeoHookeanMaterial,
	                     epsilon::SMatrix{2,2,Float64}, F, J, ip)
   sigma    .= (1.0/J)*( mat.mu*(F*F'-Identity) + mat.lambda*log(J)*Identity )
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
	k1     ::Float64         # hardening modulus

	pstrain::Vector{SVector{3,Float64}}  # plastic strain
	alpha  ::Vector{Float64}             # equivalent plastic strain
	vmStr  ::Vector{Float64}             # equivalent von Mises stress

	epsilon_dev::MVector{3,Float64}
	s_trial    ::MVector{3,Float64}
	sigma_trial::MVector{3,Float64}
	dev_Sig    ::MVector{3,Float64}

	function ElastoPlasticMaterial(E,nu,density,fy,k1,parCount)

		lambda = E*nu/(1+nu)/(1-2*nu)
		mu     = E/2/(1+nu)
		kappa  = lambda + mu

		alpha    = fill(0,parCount)
		vmStr    = fill(0,parCount)
		pstrain  = fill(zeros(3),parCount)

        new(E,nu,density,lambda,mu,kappa,fy,k1,pstrain,alpha,vmStr,
		zeros(3),zeros(3),zeros(3),zeros(3))
    end
end

function update_stress!(sig::MMatrix{2,2,Float64},mat::ElastoPlasticMaterial,
	                    eps0::SMatrix{2,2,Float64},F, J, ip)

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

export MaterialType, RigidMaterial, ElasticMaterial, ElastoPlasticMaterial, NeoHookeanMaterial
export update_stress!, computeCrackDrivingForce

end
