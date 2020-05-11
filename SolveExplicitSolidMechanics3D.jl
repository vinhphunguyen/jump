function solve_explicit_dynamics_3D(grid,solids,basis,alg::MUSL,output,Tf,dtime)
	t        = 0.
    counter  = 0
    Identity = SMatrix{2,2}(1, 0, 0, 0, 1, 0, 0, 0, 1)
	grid          = problem.grid
	solids        = problem.solids
	basis         = problem.basis
	solidCount    = problem.solidCount
	nodalMass     = grid.mass
	nodalMomentum  = grid.momentum
	nodalMomentum2 = grid.momentum2
	nodalForce    = grid.force
	interval      = problem.output.interval

	# pre_allocating arrays for temporary variables

	vel_grad      = zeros(Float64,2,2)
	D             = zeros(Float64,2,2)

	nearPoints,funcs, ders = initialise(basis,grid)

    if ( typeof(problem.output) <: PyPlotOutput )
	  pyFig_RealTime = PyPlot.figure(problem.output.figTitle,
                               figsize=problem.output.figSize, edgecolor="white", facecolor="white")
    end

	while t < Tf
	    #@printf(“Solving step: %f \n”, t)
	    # ===========================================
	    # reset grid data
	    # ===========================================
	    @inbounds for i = 1:grid.nodeCount
		  @inbounds nodalMass[i]      = 0.
		  @inbounds nodalMomentum[i]  = [0. 0. 0.]
		  @inbounds nodalMomentum2[i] = [0. 0. 0.]
		  @inbounds nodalForce[i]     = [0. 0. 0.]
	    end
	    # ===========================================
	    # particle to grid
	    # ===========================================
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
				@inbounds for i = 1:support
					in    = nearPoints[i]; # index of node ‘i’
					Ni    = funcs[i]
					dNi   = ders[:,i]
					# mass, momentum, internal force and external force
					nodalMass[in]      += Ni * fMass
					nodalMomentum[in]  += Ni * fMass * vp
					nodalForce[in]     -=      fVolume * sigma * dNi
				end
		  	end
		end

		# ===========================================
		# update grid
		# ===========================================
		@inbounds for i=1:grid.nodeCount
			nodalMomentum[i] += nodalForce[i] * dtime
	        # apply Dirichet boundary conditions
	        if grid.fixedXNodes[i] == 1
	        	nodalMomentum[i][1]  = 0.
	        	nodalForce[i][1] = 0.
	        end
	        if grid.fixedYNodes[i] == 1
	        	nodalMomentum[i][2]  = 0.
	        	nodalForce[i][2] = 0.
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
			        getShapeAndGradient(nearPoints,funcs,ders,xx[ip], grid)
					for i in 1:length(nearPoints)
						in = nearPoints[i]; # index of node ‘i’
						Ni = funcs[i]
						dNi= ders[:,i]
						vv[ip]   += (Ni * nodalForce[in] / nodalMass[in]) * dtime
					end
					# mapping the updated particle vel back to the node
					for i in 1:length(nearPoints)
						in = nearPoints[i]; # index of node ‘i’
						nodalMomentum2[in]  += funcs[i] * mm[ip] * vv[ip]
					end
			  	end
	    end
	    # # apply Dirichet boundary conditions
	    @inbounds for i = 1:grid.nodeCount
	        if grid.fixedXNodes[i] == 1
	        	nodalMomentum2[i][1]  = 0.
	        end
	        if grid.fixedYNodes[i] == 1
	        	nodalMomentum2[i][2]  = 0.
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
				support   = getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid, basis)
		        vel_grad = zeros(Float64,2,2)
				for i =1:support
					in = nearPoints[i]; # index of node ‘i’
					Ni = funcs[i]
					dNi= ders[:,i]
					vI        = nodalMomentum2[in] / nodalMass[in]
					xx[ip]   += (Ni * nodalMomentum[in] / nodalMass[in]) * dtime
					vel_grad += vI*dNi'
				end
	            D           = 0.5 * (vel_grad+vel_grad')
	            strain[ip] += dtime * D
				F[ip]      *= (Identity + vel_grad*dtime)
				vol[ip]     = det(F[ip]) * vol0[ip]
				update_stress!(stress[ip],mat,strain[ip],ip)
				#stress[ip] .= update_stress(mat,strain[ip]) #mat.lambda * (strain[ip][1,1]+strain[ip][2,2]) * Identity + 2.0 * mat.mu * strain[ip]
		  	end
		end

		if (counter%interval == 0)
			plotParticles(output,solids,[grid.lx, grid.ly],[grid.nodeCountX, grid.nodeCountY],counter)
			se,ke    = computeEnergies(solids)
			push!(problem.kinEnergy,ke)
			push!(problem.strEnergy,se)
			push!(problem.recordTime,t)
	    end

        t       += dtime
        counter += 1
    end
end

function solve_explicit_dynamics_3D(grid,solids,basis,alg::USL,output,dtime)
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

	nearPoints,funcs, ders = initialise(basis,grid)

	pyFig_RealTime = PyPlot.figure("MPM 2Disk Real-time",
                               figsize=(8/2.54, 8/2.54), edgecolor="white", facecolor="white")

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
			update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
			#stress[ip] .= update_stress(mat,strain[ip]) #mat.lambda * (strain[ip][1,1]+strain[ip][2,2]) * Identity + 2.0 * mat.mu * strain[ip]
	  	end
	end


	if (counter%interval == 0)
		plotParticles(output,solids,[grid.lx, grid.ly, grid.lz],[grid.nodeCountX, grid.nodeCountY, grid.nodeCount],counter)
		se,ke    = computeEnergies(solids)
		push!(problem.kinEnergy,ke)
		push!(problem.strEnergy,se)
		push!(problem.recordTime,t)
	end

    t       += dtime
    counter += 1
  end # end of time loop
  #@printf("Solving done \n")
end # end solve()
