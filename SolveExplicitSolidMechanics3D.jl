function solve_explicit_dynamics_3D(grid,solids,basis,alg::MUSL,output,fixes,Tf,dtime)
	t        = 0.
    counter  = 0
    Identity = UniformScaling(1.)

	solidCount    = length(solids)
	alpha         = alg.alpha

	nodalMass      = grid.mass
	nodalMomentum0 = grid.momentum0
	nodalMomentum  = grid.momentum
	nodalMomentum2 = grid.momentum2
	nodalForce     = grid.force

	# pre_allocating arrays for temporary variables


	D        = SMatrix{3,3}(0.,0.,0.,0.,0.,0.,0.,0.,0.) #zeros(Float64,2,2)

	nearPoints,funcs, ders = initialise(grid,basis)

	while t < Tf
	    #@printf(“Solving step: %f \n”, t)
	    # ===========================================
	    # reset grid data
	    # ===========================================
	    @inbounds for i = 1:grid.nodeCount
		  @inbounds nodalMass[i]      = 0.
		  @inbounds nodalMomentum0[i]  = @SVector [0., 0., 0.]
		  @inbounds nodalMomentum2[i] = @SVector [0., 0., 0.]
		  @inbounds nodalForce[i]     = @SVector [0., 0., 0.]
	    end
	    # ===========================================
	    # particle to grid
	    # ===========================================
		for s = 1:solidCount
			solid  = solids[s]
			#xx     = solid.pos
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
		        #println(nearPoints)
				@inbounds for i = 1:support
					id    = nearPoints[i]; # index of node ‘i’
					Ni    = funcs[i]
					dNix  = ders[1,i]
					dNiy  = ders[2,i]
					dNiz  = ders[3,i]
					Nim   = Ni * fMass
					# mass, momentum, internal force and external force
					nodalMass[id]      += Nim
					nodalMomentum0[id] += Nim * vp
					nodalForce[id]     -= fVolume * @SVector[sigma[1,1] * dNix + sigma[1,2] * dNiy+ + sigma[1,3] * dNiz,
					                                         sigma[1,2] * dNix + sigma[2,2] * dNiy+ + sigma[2,3] * dNiz,
					                                         sigma[1,3] * dNix + sigma[2,3] * dNiy+ + sigma[3,3] * dNiz]
				end
		  	end
		end

		# ===========================================
		# update grid
		# ===========================================
		@inbounds for i=1:grid.nodeCount
			nodalMomentum[i] = nodalMomentum0[i] + nodalForce[i] * dtime
	        # apply Dirichet boundary conditions
			if grid.fixedXNodes[i] == 1
    			nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,1)
       		    nodalMomentum[i]  = setindex(nodalMomentum[i], 0.,1)
            end
            if grid.fixedYNodes[i] == 1
    			nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,2)
    			nodalMomentum[i]  = setindex(nodalMomentum[i], 0.,2)
            end
	        if grid.fixedZNodes[i] == 1
				nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,3)
				nodalMomentum[i]  = setindex(nodalMomentum[i], 0.,3)
	        end
		end
	    # ===========================================
	    # grid to particle 1: particle vel update only
	    # ===========================================
		for s = 1:solidCount
			  	solid = solids[s]
			  	mm    = solid.mass
			  	vv    = solid.velocity
			  	@inbounds for ip = 1:solid.parCount
					vp0    = vv[ip]
					mp     = mm[ip]
					vx     = 0.
					vy     = 0.
					vz     = 0.
					support = getShapeFunctions(nearPoints,funcs,ip, grid, solid,basis)
				    for i =1:support
					   	 id = nearPoints[i]; # index of node ‘i’
					   	 Ni = funcs[i]
					   	 mI = nodalMass[id]
					   	 if ( mI > 0 )
							 mii = (1/mI)*Ni
					   		 vx += (nodalMomentum[id][1] - alpha*nodalMomentum0[id][1])*mii
					   		 vy += (nodalMomentum[id][2] - alpha*nodalMomentum0[id][2])*mii
					   		 vz += (nodalMomentum[id][3] - alpha*nodalMomentum0[id][3])*mii
					   	 end
				    end
				    vv[ip] = setindex(vv[ip],alpha*vp0[1] + vx,1)
				    vv[ip] = setindex(vv[ip],alpha*vp0[2] + vy,2)
				    vv[ip] = setindex(vv[ip],alpha*vp0[3] + vz,3)
					# mapping the updated particle vel back to the node
					for i in 1:support
						id = nearPoints[i] # index of node ‘i’
						nodalMomentum2[id]  += funcs[i] * mp * vv[ip]
					end
			  	end
	    end
	    # # apply Dirichet boundary conditions
	    @inbounds for i = 1:grid.nodeCount
	        if grid.fixedXNodes[i] == 1
				nodalMomentum2[i] = setindex(nodalMomentum2[i],0.,1)
	        end
	        if grid.fixedYNodes[i] == 1
	        	nodalMomentum2[i] = setindex(nodalMomentum2[i],0.,2)
	        end
			if grid.fixedZNodes[i] == 1
			    nodalMomentum2[i] = setindex(nodalMomentum2[i],0.,3)
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
		  	stress = solid.stress
		  	strain = solid.strain
		  	mat    = solid.mat
		  	@inbounds for ip = 1:solid.parCount
				support   = getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid, basis)
				vel_grad  = SMatrix{3,3}(0., 0., 0., 0.,0,0,0,0,0)
				xxp       = xx[ip]
				for i = 1:support
					id   = nearPoints[i] # index of node ‘i’
					Ni   = funcs[i]
					dNix = ders[1,i]
					dNiy = ders[2,i]
					dNiz = ders[3,i]
					m    = nodalMass[id]
					if ( m > 0.)
						mii = 1/m
						vIx        = mii * nodalMomentum2[id][1]
						vIy        = mii * nodalMomentum2[id][2]
						vIz        = mii * nodalMomentum2[id][3]
					    xxp       += Ni  * mii * dtime * nodalMomentum[id]
				        vel_grad  += SMatrix{3,3}(dNix*vIx, dNiy*vIx, dNiz*vIx,
				                                  dNix*vIy, dNiy*vIy, dNiz*vIy,
												  dNix*vIz, dNiy*vIz, dNiz*vIz)
				    end
				end
				xx[ip]      = xxp
	            D           = 0.5 * (vel_grad + vel_grad')
	            strain[ip]  += dtime * D
				F[ip]       *= (Identity + vel_grad*dtime)
				J            = det(F[ip])
				vol[ip]      = J * vol0[ip]

				if ( J < 0. )
					@printf("Troubled particle: %f %f \n", xx[ip][1], xx[ip][2])
					println(F[ip])
					closeFile(fixes)
					@error("J is negative\n")
			    end

				update_stress!(stress[ip],mat,strain[ip],dtime*D,F[ip],J,ip,dtime)
		  	end
		end

		if (counter%output.interval == 0)
            plotParticles_3D(output,solids,[grid.lx, grid.ly, grid.lz],
			            [grid.nodeCountX, grid.nodeCountY, grid.nodeCount],counter)
			compute(fixes,t)
	    end

        t       += dtime
        counter += 1
    end
	closeFile(fixes)
end

######################################################################
#  Update Stress Last:: UNFINISHED!!!
######################################################################


function solve_explicit_dynamics_3D(grid,solids,basis,alg::USL,output,fixes,Tf,dtime)
    t       = 0.
    counter = 0

    Identity = SMatrix{3,3}(1., 0., 0., 0., 1., 0., 0., 0., 1.)

	solidCount    = length(solids)
	nodalMass     = grid.mass
	nodalMomentum = grid.momentum
	nodalForce    = grid.force

	interval      = output.interval

	# pre_allocating arrays for temporary variables

	vel_grad      = zeros(Float64,3,3)
	D             = zeros(Float64,3,3)

	nearPoints,funcs, ders = initialise(grid,basis)

  while t < Tf

    #@printf("Solving step: %f \n", t)

    # ===========================================
    # reset grid data
    # ===========================================

    @inbounds for i = 1:grid.nodeCount
	  nodalMass[i]      = 0.
	  nodalMomentum[i]  = [0. 0. 0.]
	  nodalForce[i]     = [0. 0. 0.]
    end

    # ===========================================
    # particle to grid
    # ===========================================

	@inbounds  for s = 1:solidCount
		solid  = solids[s]
		xx     = solid.pos
		mm     = solid.mass
		vv     = solid.velocity
		vol    = solid.volume
		stress = solid.stress

	  	@inbounds  for ip = 1:solid.parCount
			support   = getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid, basis)
	        fVolume   = vol[ip]
	        fMass     = mm[ip]
	        vp        = vv[ip]
	        sigma     = stress[ip]
#println(nearPoints)
			for i = 1:support
				in    = nearPoints[i]; # index of node 'i'
				Ni    = funcs[i]
				dNi   = ders[:,i]
				# mass, momentum, internal force and external force
				nodalMass[in]      += Ni * fMass
				nodalMomentum[in]  += Ni * fMass * vp
				nodalForce[in]     -=      fVolume * sigma * dNi
			end
	  	end
	end
#println("Done particle to grid.")
	# ===========================================
	# update grid
	# ===========================================

	@inbounds  for i=1:grid.nodeCount

		nodalMomentum[i] += nodalForce[i] * dtime
        # apply Dirichet boundary conditions
        if grid.fixedXNodes[i] == 1
        	nodalMomentum[i][1]  = 0.
        	nodalForce[i][1] = 0
        end
        if grid.fixedYNodes[i] == 1
        	nodalMomentum[i][2]  = 0.
        	nodalForce[i][2] = 0
        end
		if grid.fixedZNodes[i] == 1
			nodalMomentum[i][3]  = 0.
			nodalForce[i][3] = 0
		end
	end
    #println("HEHEHE")
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

	  	mat    = solid.mat

	  	@inbounds for ip = 1:solid.parCount
            support=getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid, basis)
	        vel_grad .= 0.
			# println(nearPoints)
			# println(vel_grad)
			# println(ders[:,1])
			for i  = 1:support
				in = nearPoints[i]; # index of node 'i'
				Ni = funcs[i]
				dNi= ders[:,i]
			    if(nodalMass[in] > alg.tolerance)
					vI        = nodalMomentum[in] / nodalMass[in]
					vv[ip]   += (Ni * nodalForce[in] / nodalMass[in]) * dtime
					xx[ip]   += (Ni * nodalMomentum[in] / nodalMass[in]) * dtime
					vel_grad += vI*dNi'
				end
			end
            D          .= 0.5 * (vel_grad+vel_grad')
            strain[ip] += dtime * D
			F[ip]      *= (Identity + vel_grad*dtime)
            J           =  det(F[ip])
			vol[ip]     = J * vol0[ip]
			
			update_stress!(stress[ip],mat,strain[ip],dtime*D,F[ip],J,ip,dtime)
			#stress[ip] .= update_stress(mat,strain[ip]) #mat.lambda * (strain[ip][1,1]+strain[ip][2,2]) * Identity + 2.0 * mat.mu * strain[ip]
	  	end
	end

	if (counter%output.interval == 0)
        plotParticles_3D(output,solids,[grid.lx, grid.ly, grid.lz],
		            [grid.nodeCountX, grid.nodeCountY, grid.nodeCount],counter)
		compute(fixes,t)
    end

    t       += dtime
    counter += 1
  end # end of time loop
  #@printf("Solving done \n")
end # end solve()
