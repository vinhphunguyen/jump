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

# This file contains functions for solving 2D explicit dynamics problems
# coupling with phase field fracture. One pass AM solver.
#

using Fix


######################################################################
# Modified Update Stress Last
######################################################################
function solve_explicit_dynamics_PFM_2D(grid,solids,basis,alg::MUSL,output,Tf,dtime)
	t       = 0.
    counter = 0
    Identity      = SMatrix{2,2}(1, 0, 0, 1)

	solidCount    = length(solids)

	nodalMass      = grid.mass
	nodalMomentum  = grid.momentum
	nodalMomentum0 = grid.momentum0
	nodalMomentum2 = grid.momentum2
	nodalForce     = grid.force
	nodalPFForce   = grid.pfForce
	nodalDamage0   = grid.damage0
	nodalDamage    = grid.damage

	vel_grad      = zeros(Float64,2,2)
	D             = zeros(Float64,2,2)
	vvp           = zeros(2)
	damgrad       = zeros(2)
	interval      = output.interval

	# allocate memory for grid basis and grads once
	nearPoints,funcs, ders = initialise(basis,grid)
	nearPointsLin    = [0, 0, 0, 0]
	funcsLin         = [0., 0., 0., 0.]
	linBasis         = LinearBasis()

	if ( typeof(basis) <: CPDIQ4Basis )
		nearPoints = Vector{Int64}(undef,16)
		[nearPoints[i]=0 for i=1:16]
	    funcs = zeros(16)
		ders  = zeros(2,16)
	end

    if ( typeof(output) <: PyPlotOutput )
	  pyFig_RealTime = PyPlot.figure(output.figTitle,
                               figsize=output.figSize, edgecolor="white", facecolor="white")
    end

	while t < Tf
	    #@printf(“Solving step: %f \n”, t)
	    # ===========================================
	    # reset grid data
	    # ===========================================
	     @inbounds for i = 1:grid.nodeCount
		  nodalMass[i]      = 0.
		  nodalMomentum0[i] = @SVector [0., 0.]
		  nodalMomentum2[i] = @SVector [0., 0.]
		  nodalForce[i]     = @SVector [0., 0.]
	    end
	    # ===========================================
	    # particle to grid
	    # ===========================================
		#println(problem.bodyforce.g)
		for s = 1:solidCount
			solid  = solids[s]
			# deformable solids only
			if solid.rigid continue end
			xx     = solid.pos
			mm     = solid.mass
			vv     = solid.velocity
			vol    = solid.volume
			stress = solid.stress
			dam    = solid.damage
		  	@inbounds for ip = 1:solid.parCount
		        support   = getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid, basis)
		        fVolume   = vol[ip]
		        fMass     = mm[ip]
		        vp        = vv[ip]
		        sigma     = stress[ip]
				d         = dam[ip]
				degrad    = (1-d)^2 # AT1/2 only
				#body      = problem.bodyforce(xx[ip],t)
				@inbounds for i = 1:support
					in    = nearPoints[i]; # index of node ‘i’
					Ni    = funcs[i]
					dNi   = ders[:,i]
					# mass, momentum, internal force and external force
					nodalMass[in]       += Ni * fMass
					nodalMomentum0[in]  += Ni * fMass * vp
					nodalForce[in]      -=      fVolume * sigma * dNi * degrad
					nodalForce[in]      +=      fMass   * body  *  Ni
				end
		  	end
		end

		# ===========================================
		# update grid
		# ===========================================
		@inbounds for i=1:grid.nodeCount
			nodalMomentum[i] = nodalMomentum0[i] + nodalForce[i] * dtime
	        # # apply Dirichet boundary conditions
	        # if grid.fixedXNodes[i] == 1
	        # 	#nodalMomentum0[i][1]  = 0.
	        # 	#nodalMomentum[i][1]  = 0.
			# 	nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,1)
			# 	nodalMomentum[i]  = setindex(nodalMomentum[i],0.,1)
	        # end
	        # if grid.fixedYNodes[i] == 1
	        # 	# nodalMomentum0[i][2]  = 0.
	        # 	# nodalMomentum[i][2]  = 0.
			#
			# 	nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,2)
			# 	nodalMomentum[i]  = setindex(nodalMomentum[i],0.,2)
	        # end
		end

		@inbounds for s = 1:solidCount
			solid  = solids[s]
			# deformable solids only
			if !solid.rigid continue end
			xx     = solid.pos
			ve     = solid.mat.vel
			@inbounds for ip = 1:solid.parCount
				getAdjacentGridPoints(nearPointsLin,xx[ip],grid,linBasis)
	#			println(nearPoints)
				@inbounds for i = 1:length(nearPointsLin)
					in                  = nearPointsLin[i]; # index of node 'i'
					nodalMomentum0[in]  = nodalMass[in] * ve
					nodalMomentum[in]   = nodalMass[in] * ve
					#println(nodalMomentum)
				end
			end
		end

	    # ===========================================
	    # grid to particle 1: particle vel update only
	    # ===========================================
		for s = 1:solidCount
			  	solid = solids[s]
				if solid.rigid continue end
			  	xx    = solid.pos
			  	mm    = solid.mass
			  	vv    = solid.velocity
			  	@inbounds @simd for ip = 1:solid.parCount
					vp0     = vv[ip]
					vx     = 0.
					vy     = 0.
			        support = getShapeFunctions(nearPoints,funcs,ip, grid, solid, basis)
					for i =1:support
						in = nearPoints[i]; # index of node ‘i’
						Ni = funcs[i]
						mI = nodalMass[in]
						if ( mI > 0 )
							vx += Ni * (nodalMomentum[in][1] - alg.alpha*nodalMomentum0[in][1])/mI
							vy += Ni * (nodalMomentum[in][2] - alg.alpha*nodalMomentum0[in][2])/mI
						end
					end
					vv[ip] = setindex(vv[ip],alg.alpha*vp0[1] + vx,1)
					vv[ip] = setindex(vv[ip],alg.alpha*vp0[2] + vy,2)
					# mapping the updated particle vel back to the node
					for i in 1:support
						in = nearPoints[i] # index of node ‘i’
						nodalMomentum2[in]  += funcs[i] * mm[ip] * vv[ip]
					end
			  	end
	    end

		@inbounds for s = 1:solidCount
   			solid  = solids[s]
   			# deformable solids only
   			if !solid.rigid continue end
   			xx     = solid.pos
   			ve     = solid.mat.vel
   			@inbounds for ip = 1:solid.parCount
   				getAdjacentGridPoints(nearPointsLin,xx[ip],grid,linBasis)
   	#			println(nearPoints)
   				@inbounds for i = 1:length(nearPointsLin)
   					in                  = nearPointsLin[i]; # index of node 'i'
   					nodalMomentum2[in]  = nodalMass[in] * ve
   					#println(nodalMomentum)
   				end
   			end
   		end

	    # ===========================================
	    # particle to grid 2:
	    # ===========================================
		@inbounds for s = 1:solidCount
		  	solid = solids[s]
			if solid.rigid continue end
		  	xx    = solid.pos
		  	mm    = solid.mass
		  	vv    = solid.velocity
		  	vol   = solid.volume
		  	vol0  = solid.volumeInitial
		  	F     = solid.deformationGradient
		  	mat   = solid.mat
		  	stress = solid.stress
		  	strain = solid.strain
			crackForce    = solid.crackForce
		  	@inbounds @simd for ip = 1:solid.parCount
				support=getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid, basis)
				#getShapeAndGradient(nearPoints,funcs,ders,xx[ip], grid)
		        vel_grad .= 0. #zeros(Float64,2,2)
				for i = 1:support
					in = nearPoints[i]; # index of node ‘i’
					Ni = funcs[i]
					dNi= ders[:,i]
					m = nodalMass[in]
					if ( m > 0.)
						vI         = nodalMomentum2[in] /m
						xx[ip]    += (Ni * nodalMomentum[in]/m) * dtime
						vel_grad .+= vI*dNi'
				    end
				end
	            D            .= 0.5 * (vel_grad.+vel_grad')
	            strain[ip]  += dtime * D
				F[ip]       *= (Identity + vel_grad*dtime)
				J           = det(F[ip])
				if ( J < 0. )
					@printf("Troubled particle: %f %f \n", xx[ip][1], xx[ip][2])
					println(F[ip])
					@error("J is negative\n")
				end
				vol[ip]     = J * vol0[ip]
				update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
				crackForce[ip] = computeCrackDrivingForce(stress[ip],mat)
		  	end
		end

		# if CPDI, do extra thing here

		if ( typeof(basis) <: CPDIQ4Basis )
			@inbounds for s = 1:solidCount
				solid = solids[s]
				@inbounds for c = 1:size(solid.nodes,2)
				  getShapeFuncs(nearPointsLin,funcsLin, solid.nodes[:,c], grid, solid, LinearBasis())
				  for i = 1:length(nearPointsLin)
					  in = nearPointsLin[i]; # index of node ‘i’
					  Ni = funcsLin[i]
					  #vI        = nodalMomentum2[in] / nodalMass[in]
					  if nodalMass[in] > 0.
						  solid.nodes[:,c]   += (Ni * nodalMomentum[in] / nodalMass[in]) * dtime
					  end
				  end
			    end
		    end
	    end

		##################################################
		# update position of rigid solids
		##################################################
		@inbounds for s = 1:solidCount
			# only rigid solids here
			solid = solids[s]
			if !solid.rigid continue end

			xx    = solid.pos
			ve    = solid.mat.vel
			#println(ve)
			@inbounds for ip = 1:solid.parCount
			  xx[ip]   += ve * dtime
			end
		end

		############################################################
		# solving the damage sub-problem with updated crack driving
		# force
		############################################################
		@inbounds for i = 1:grid.nodeCount
		 nodalMass[i]      = 0.
		 nodalDamage0[i]   = 0.
		 nodalPFForce[i]   = 0.
	   end
        xi = dtime/0.00057000
		for s = 1:solidCount
			solid  = solids[s]
			if solid.rigid continue end
			xx     = solid.pos
			mm     = solid.mass
			vol    = solid.volume
			dam    = solid.damage
			crackForce    = solid.crackForce
			gradDamage= solid.gradDamage
			Gf     = solid.mat.Gf
			b      = solid.mat.l0
			calpha = 2.666666667
			@inbounds for ip = 1:solid.parCount
				support   = getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid, basis)
				fVolume   = vol[ip]
				d         = dam[ip]
				Ybar      = crackForce[ip]
				gradDam   = gradDamage[ip]
				dalpha    = 1.0 # AT1 model
				domega    = 2.0 * ( d - 1.0)
				coef1     = Gf/calpha/b
				coef2     = 2*b*Gf/calpha
				@inbounds for i = 1:support
					in    = nearPoints[i]; # index of node ‘i’
					Ni    = funcs[i]
					dNi   = ders[:,i]
					# mass, momentum, internal force and external force
					nodalMass[in]     += Ni  * fVolume
					nodalDamage0[in]  += Ni  * fVolume * d
					nodalPFForce[in]  += (Ni * (domega*Ybar + dalpha*coef1) +
					                        coef2*(dNi[1]*gradDam[1]+dNi[2]*gradDam[2])) * fVolume
				end
			end
		end

		# ===========================================
		# update grid
		# ===========================================
		@inbounds for i=1:grid.nodeCount
			dI = nodalDamage0[i] - xi*nodalPFForce[i]
			if dI < 0   dI = 0. end
			if dI > 1.0 dI = 1. end
			nodalDamage[i] = dI
			# apply Dirichet boundary conditions
			# if grid.fixedXNodes[i] == 1
			# 	#nodalMomentum0[i][1]  = 0.
			# 	#nodalMomentum[i][1]  = 0.
			# 	nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,1)
			# 	nodalMomentum[i]  = setindex(nodalMomentum[i],0.,1)
			# end
			# if grid.fixedYNodes[i] == 1
			# 	# nodalMomentum0[i][2]  = 0.
			# 	# nodalMomentum[i][2]  = 0.
			#
			# 	nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,2)
			# 	nodalMomentum[i]  = setindex(nodalMomentum[i],0.,2)
			# end
		end

		# ===========================================
	    # particle to grid: update particle damage
	    # ===========================================
		@inbounds for s = 1:solidCount
		  	solid = solids[s]
		  	xx    = solid.pos
		  	mm    = solid.mass
		  	vol   = solid.volume
			dam   = solid.damage
			gradDamage= solid.gradDamage
		  	@inbounds @simd for ip = 1:solid.parCount
				support  = getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid, basis)
				damgrad .= 0.
				for i = 1:support
					in = nearPoints[i]; # index of node ‘i’
					Ni = funcs[i]
					dNi= ders[:,i]
					m =  nodalMass[in]
					if ( m > 0.)
						dam[ip]  += Ni * ( nodalDamage[in] - nodalDamage0[in]) / m
						damgrad  += dNi*nodalDamage[in] /m
				    end
				end
				gradDamage[ip] = damgrad
		  	end
		end

		if (counter%interval == 0)
			plotParticles(output,solids,[grid.lx, grid.ly],
			             [grid.nodeCountX, grid.nodeCountY],counter)
			se,ke    = computeEnergies(solids)
			# push!(problem.kinEnergy,ke)
			# push!(problem.strEnergy,se)
			# push!(problem.recordTime,t)
			# [compute(fix,t) for fix in problem.fixes]
	    end

        t       += dtime
        counter += 1
    end
	#[closeFile(fix) for fix in problem.fixes]
end
