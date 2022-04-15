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

# This file contains functions for solving 1D explicit dynamics problems
# USL and MUSL are provided, each function for each algorithm
#

using Fix

function solve_explicit_dynamics_1D(grid,solids,basis,alg::MUSL,interval,Tf,dtime,bodyforce)
	t       = 0.
  counter = 0

	solidCount     = length(solids)
	nodalMass      = grid.mass
	nodalMomentum  = grid.momentum
	nodalMomentum2 = grid.momentum2
	nodalForce     = grid.force
	#interval      = problem.output.interval

	# allocate memory for grid basis and grads once
	nearPoints, funcs, ders = initialise(grid, basis)


	while t < Tf
        t       += dtime
        counter += 1

	    #@printf(“Solving step: %f \n”, t)
	    # ===========================================
	    # reset grid data
	    # ===========================================
	     @inbounds for i = 1:grid.nodeCount
		  nodalMass[i]      = 0.
		  nodalMomentum[i]  = 0.
		  nodalMomentum2[i] = 0.
		  nodalForce[i]     = 0.
	    end
	    # ===========================================
	    # particle to grid
	    # ===========================================
		#println(problem.bodyforce.g)
		for s = 1:solidCount
			solid  = solids[s]
			xx     = solid.pos
			mm     = solid.mass
			vv     = solid.velocity
			vol    = solid.volume
			stress = solid.stress
		  	@inbounds for ip = 1:solid.parCount
		        support   = getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid, basis)
		        fVolume   = vol[ip]
		        fMass     = mm[ip]
		        vp        = vv[ip]
		        sigma     = stress[ip]
				    body      = bodyforce(xx[ip],t)
				@inbounds for i = 1:support
					in    = nearPoints[i]; # index of node ‘i’
					Ni    = funcs[i]
					dNi   = ders[i]
					# mass, momentum, internal force and external force
					nodalMass[in]      += Ni * fMass
					nodalMomentum[in]  += Ni * fMass * vp
					nodalForce[in]     -=      fVolume * sigma * dNi
					nodalForce[in]     +=      fMass   * body  *  Ni
				end
		  	end
		end

		# ===========================================
		# update grid
		# ===========================================
		@inbounds for i=1:grid.nodeCount
			nodalMomentum[i] += nodalForce[i] * dtime
	        # apply Dirichet boundary conditions
	        if grid.fixedNodes[i] == 1
	        	nodalMomentum[i]  = 0.
	        	nodalForce[i]     = 0.
	        end
		end
	    # ===========================================
	    # grid to particle 1: particle vel update only
	    # ===========================================
		for s = 1:solidCount
			  	solid = solids[s]
			  	xx    = solid.pos
			  	mm    = solid.mass
			  	vv    = solid.velocity
			  	vol   = solid.volume
			  	vol0  = solid.volumeInitial
			  	F     = solid.deformationGradient
			  	mat   = solid.mat
			  	stress = solid.stress
			  	strain = solid.strain
			  	mat    = solid.mat
			  	@inbounds for ip = 1:solid.parCount
			        support = getShapeFunctions(nearPoints,funcs,ip, grid, solid, basis)
					for i =1:support
						in = nearPoints[i]; # index of node ‘i’
						Ni = funcs[i]
						if ( nodalMass[in] > 0 ) vv[ip]   += (Ni * nodalForce[in] / nodalMass[in]) * dtime end
					end
					# mapping the updated particle vel back to the node
					for i in 1:support
						in = nearPoints[i]; # index of node ‘i’
						nodalMomentum2[in]  += funcs[i] * mm[ip] * vv[ip]
					end
			  	end
	    end
	    # # apply Dirichet boundary conditions
	    @inbounds for i = 1:grid.nodeCount
	        if grid.fixedNodes[i] == 1
	        	nodalMomentum2[i]  = 0.
	        end
	    end

	    # ===========================================
	    # particle to grid 2:
	    # ===========================================
		for s = 1:solidCount
		  	solid = solids[s]
		  	xx    = solid.pos
		  	mm    = solid.mass
		  	vv    = solid.velocity
		  	vol   = solid.volume
		  	vol0  = solid.volumeInitial
		  	F     = solid.deformationGradient
		  	mat   = solid.mat
		  	stress = solid.stress
		  	strain = solid.strain
		  	mat    = solid.mat
		  	@inbounds for ip = 1:solid.parCount
				support=getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid, basis)
		        vel_grad = 0.
				for i = 1:support
					in = nearPoints[i]; # index of node ‘i’
					Ni = funcs[i]
					dNi= ders[i]
					m = nodalMass[in]
					if ( m > 0.)
						vI        = nodalMomentum2[in] /m
						xx[ip]   += (Ni * nodalMomentum[in]/m) * dtime
						vel_grad += vI*dNi
				    end
				end
	            D           = vel_grad
	            strain[ip] += dtime * D
				F[ip]      *= (1 + vel_grad*dtime)
				J           = F[ip]
				if ( J < 0. )
					@printf("Troubled particle: %f %f \n", xx[ip][1], xx[ip][2])
					println(F[ip])
					@error("J is negative\n")
				end
				vol[ip]     = J * vol0[ip]
				stress[ip]  = mat.lambda*log(J)/J + mat.mu*J - mat.mu/J#; % Neo-Hookean
				#update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
				#stress[ip] .= update_stress(mat,strain[ip]) #mat.lambda * (strain[ip][1,1]+strain[ip][2,2]) * Identity + 2.0 * mat.mu * strain[ip]
		  	end
		end

		if (counter%interval == 0)
			# push!(problem.recordTime,t)
            # problem.recordX = [problem.recordX solids[1].pos]
		end

    end
end

function solve_explicit_dynamics_1D(grid,solids,basis,alg::USL,interval,Tf,dtime,bodyforce)
    t       = 0.
    counter = 0

	solidCount    = length(solids)
	nodalMass     = grid.mass
	nodalMomentum = grid.momentum
	nodalForce    = grid.force

	nearPoints, funcs, ders = initialise(grid,basis)


  while t < Tf

    #@printf("Solving step: %f \n", t)

    # ===========================================
    # reset grid data
    # ===========================================

    @inbounds for i = 1:grid.nodeCount
	  nodalMass[i]      = 0.
	  nodalMomentum[i]  = 0.
	  nodalForce[i]     = 0.
    end

    # ===========================================
    # particle to grid
    # ===========================================

	@inbounds for s = 1:solidCount
		solid  = solids[s]
		xx     = solid.pos
		mm     = solid.mass
		vv     = solid.velocity
		vol    = solid.volume
		stress = solid.stress
		#gradVe = solid.gradVelo

	  	@inbounds for ip = 1:solid.parCount
	        #getShapeAndGradient(nearPoints,funcs,ders,xx[ip], grid)
			support=getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid, basis)

	        fVolume   = vol[ip]
	        fMass     = mm[ip]
	        vp        = vv[ip]
	        sigma     = stress[ip]
			#vgrad     = gradVe[ip]
			    body      = bodyforce(xx[ip],t)
			#println(nearPoints)
			@inbounds for i = 1:support
				in    = nearPoints[i]; # index of node 'i'
				Ni    = funcs[i]
				dNi   = ders[i]
				# mass, momentum, internal force and external force
				nodalMass[in]      += Ni * fMass
				nodalMomentum[in]  += Ni * fMass * (vp)# + vgrad*(grid.pos[in]-xx[ip]))
				nodalForce[in]     -=      fVolume * sigma * dNi
				nodalForce[in]     +=      fMass   * body  *  Ni
			end
	  	end
	end

	# ===========================================
	# update grid
	# ===========================================

	@inbounds for i=1:grid.nodeCount

		nodalMomentum[i] += nodalForce[i] * dtime
        # apply Dirichet boundary conditions
        if grid.fixedNodes[i] == 1
        	nodalMomentum[i] = 0.
        	nodalForce[i]    = 0
        end
	end

    # ===========================================
    # grid to particle
    # ===========================================

    @inbounds for s = 1:solidCount
	  	solid = solids[s]
	  	xx    = solid.pos
	  	mm    = solid.mass
	  	vv    = solid.velocity
	  	vol   = solid.volume
	  	vol0  = solid.volumeInitial
	  	F     = solid.deformationGradient
	  	mat   = solid.mat
	  	stress = solid.stress
	  	strain = solid.strain
		#vgrad  = solid.gradVelo

	  	mat    = solid.mat

	  	@inbounds for ip = 1:solid.parCount
			support  = getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid, basis)
	        vel_grad = 0
			for i = 1:support
				in = nearPoints[i]; # index of node 'i'
				Ni = funcs[i]
				dNi= ders[i]
			    if(nodalMass[in] > alg.tolerance)
					vI        = nodalMomentum[in] / nodalMass[in]
					vv[ip]   += (Ni * nodalForce[in] / nodalMass[in]) * dtime
					xx[ip]   += (Ni * nodalMomentum[in] / nodalMass[in]) * dtime
					vel_grad += vI*dNi
				end
			end
            D           = vel_grad
            strain[ip] += dtime * D
			F[ip]      *= (1 + vel_grad*dtime)
			J           = F[ip]
			vol[ip]     = J * vol0[ip]
            stress[ip]  = mat.lambda*log(J)/J + mat.mu*J - mat.mu/J#; % Neo-Hookean
			#update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
			#stress[ip] .= update_stress(mat,strain[ip]) #mat.lambda * (strain[ip][1,1]+strain[ip][2,2]) * Identity + 2.0 * mat.mu * strain[ip]
	  	end
	end


    t       += dtime
    counter += 1

	if (counter%interval == 0)
		# if (counter%interval == 0)
		# 	push!(problem.recordTime,t)
        #     problem.recordX = [problem.recordX solids[1].pos]
		# end
	end

  end # end of time loop
  #@printf("Solving done \n")
  #[closeFile(fix) for fix in problem.fixes]
end # end solve()
