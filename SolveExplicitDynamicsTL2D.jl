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
# using TLMPM. USL and MUSL are provided, each function for each algorithm
#

using Fix


######################################################################
# Update Stress Last
######################################################################
function solve_explicit_dynamics_2D_TL(grid,solids,basis,alg::USL,output,Tf,dtime)
    t       = 0.
    counter = 0

    Identity = SMatrix{2,2}(1, 0, 0, 1)


	solidCount     = length(solids)
    interval       = output.interval

	linBasis       = LinearBasis()
	nearPointsLin  = [0,0,0,0]

    # pre_allocating arrays for temporary variables

 	Fdot          = zeros(Float64,2,2)
	D             = zeros(Float64,2,2)

	nearPoints,funcs, ders = initialise(basis,grids[1])

	if ( typeof(basis) <: CPDIQ4Basis )
		nearPointsLin    = [0, 0, 0, 0]
		funcsLin         = [0., 0., 0., 0.]
	end


    if ( typeof(output) <: PyPlotOutput )
	  pyFig_RealTime = PyPlot.figure(output.figTitle,
                               figsize=output.figSize, edgecolor="white", facecolor="white")
    end

    # compute once nodal mass
	# and probably shape functions ... and stored
	# and color field to recognise surface particles
	@inbounds for s = 1:solidCount
  	  solid  = solids[s]
	  grid   = grids[s]
	  # detect surface particles
	  mpoints = doCellParticleInteraction(solid,grid)
  	  # deformable solids only
  	  if solid.rigid continue end
  	  XX        = solid.pos0
  	  mm        = solid.mass
	  #rad       = solid.radius
	  vol       = solid.volumeInitial
  	  color     = solid.color
	  nodalMass = grid.mass
	  rp        = 1.0 / sqrt(grid.dx^2+grid.dy^2)
  	  @inbounds for ip = 1:solid.parCount
  		  support   = getShapeAndGradientTL(nearPoints,funcs,ders,ip, grid, solid, basis)
  		  fMass     = mm[ip]
		  #rad[ip]   = 0.5*cbrt(vol[ip])  # radius of the sphere of particle 'ip'
  		  #println(nearPoints)
  		  @inbounds for i = 1:support
  			  in              = nearPoints[i]; # index of node 'i'
  			  nodalMass[in]  += funcs[i] * fMass
  		  end
		  ################################################
		  # find surface particles by computing the sum kernel
		  xp = XX[ip][1]
		  yp = XX[ip][2]
		  e  = floor(Int64,(xp-grid.xmin)/grid.dx) + 1 +
		                   (grid.nodeCountX-1)*floor(Int64,(yp-grid.ymin)/grid.dy)
		  neighbors=getNeighbors(e, grid)
		  for ie=1:length(neighbors)
			  elemId = neighbors[ie]
			  pars   = mpoints[elemId]
			  #println(pars)
			  for iq = 1:length(pars)
				  q  = pars[iq]
				  xq = XX[q][1]
				  yq = XX[q][2]
				  dist = sqrt( (xp-xq)^2 +(yp-yq)^2 ) * rp
				  #println(dist)
                  if dist <= 1.
					  kernel     = 1-3*dist^2+2*dist^3
					  color[ip] += kernel
				  end
			  end
		  end
  	  end
    end

  while t < Tf

    @printf("Solving step: %d %f \n", counter, t)

    # ===========================================
    # reset grid data
    # ===========================================

	@inbounds for g=1:solidCount
		grid           = grids[g]
		nodalMomentum0 = grid.momentum0
		nodalForce     = grid.force
	    @inbounds for i = 1:grid.nodeCount
		  nodalMomentum0[i] = [0. 0.]
		  nodalForce[i]     = [0. 0.]
	    end
    end

    # ===========================================
    # particle to grid (deformable solids)
    # ===========================================

	@inbounds for s = 1:solidCount
		solid  = solids[s]
		grid   = grids[s]
		# deformable solids only
		if solid.rigid continue end

		nodalMomentum0 = grid.momentum0
		nodalForce     = grid.force

		XX     = solid.pos0
		mm     = solid.mass
		vv     = solid.velocity
		vol    = solid.volume
		stress = solid.stress
		F      = solid.deformationGradient
		vol0   = solid.volumeInitial

	  	@inbounds for ip = 1:solid.parCount
			support   = getShapeAndGradientTL(nearPoints,funcs,ders,ip, grid, solid, basis)
	        fVolume   = vol0[ip]
	        fMass     = mm[ip]
	        vp        = vv[ip]
	        sigma     = stress[ip]
			P         = det(F[ip])*sigma*inv(F[ip])'
			#body      = problem.bodyforce(XX[ip],t)
			#println(nearPoints)
			@inbounds for i = 1:support
				in    = nearPoints[i]; # index of node 'i'
				Ni    = funcs[i]
				dNi   = ders[:,i]
				# momentum, internal force and external force and contact force
				nodalMomentum0[in] += Ni * fMass * vp #+ vgrad*(grid.pos[in]-xx[ip]))
				nodalForce[in]     -=      fVolume * P     * dNi
				nodalForce[in]     +=      fMass   * body  *  Ni
			end
	  	end
	end

	# ===========================================
	# update grids (each grid for each solid)
	# ===========================================

	@inbounds for g=1:solidCount
		grid           = grids[g]
		nodalMomentum0 = grid.momentum0
		nodalMomentum  = grid.momentum
		nodalForce     = grid.force
		@inbounds for i=1:grid.nodeCount
			nodalMomentum[i] .= nodalMomentum0[i] + nodalForce[i] * dtime
	        # apply Dirichet boundary conditions
	        if grid.fixedXNodes[i] == 1
	        	nodalMomentum0[i][1]  = 0.
	        	nodalMomentum[i][1]   = 0.
	        end
	        if grid.fixedYNodes[i] == 1
	        	nodalMomentum0[i][2] = 0.
	        	nodalMomentum[i][2]  = 0.
	        end
		end
    end

    # ===========================================
    # grid to particle (deformable solids)
    # ===========================================

    @inbounds for s = 1:solidCount
	  	solid = solids[s]
		grid  = grids[s]

		nodalMass      = grid.mass
		nodalMomentum  = grid.momentum
		nodalMomentum0  = grid.momentum0

	  	xx    = solid.pos
	  	mm    = solid.mass
	  	vv    = solid.velocity
	  	vol   = solid.volume
	  	vol0  = solid.volumeInitial
	  	F     = solid.deformationGradient
	  	mat   = solid.mat
	  	stress = solid.stress
	  	strain = solid.strain
		vgrad  = solid.gradVelo

	  	@inbounds for ip = 1:solid.parCount
			support   = getShapeAndGradientTL(nearPoints,funcs,ders,ip, grid, solid, basis)
	        Fdot     .= 0.
			for i = 1:support
				in = nearPoints[i]; # index of node 'i'
				Ni = funcs[i]
				dNi= ders[:,i]
				mI = nodalMass[in]
			    if (mI > alg.tolerance)
					invM      = 1.0 / mI
					vI        = nodalMomentum[in] * invM
					vv[ip]   += Ni * (nodalMomentum[in] - nodalMomentum0[in]) * invM
					xx[ip]   += Ni * vI * dtime
					Fdot     += vI*dNi'
				end
			end
			F[ip]      += Fdot*dtime
			vel_grad    = Fdot*inv(F[ip])
			D           .= 0.5 * (vel_grad+vel_grad')
			strain[ip] += dtime * D
			J           = det(F[ip])
			vol[ip]     = J * vol0[ip]
			update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
	  	end
	end

    #### CPDI specific #################################
	### update corners
	if ( typeof(basis) <: CPDIQ4Basis )
		@inbounds for s = 1:solidCount
			solid = solids[s]
			if solid.rigid continue end
			@inbounds for c = 1:size(solid.nodes,1)
			  getShapeFuncs(nearPointsLin,funcsLin, solid.nodes[:,c], grid, solid, linBasis)
			  for i in 1:length(nearPointsLin)
				  in = nearPointsLin[i]; # index of node ‘i’
				  Ni = funcsLin[i]
				  #vI        = nodalMomentum2[in] / nodalMass[in]
				  if nodalMass[in] > alg.tolerance
					  solid.nodes[:,c]  += (Ni * nodalMomentum[in] / nodalMass[in]) * dtime
				  end
			  end
			end
		end
	end


	if (counter%interval == 0)
		plotParticles(output,solids,[grids[1].lx, grids[1].ly],
		              [grids[1].nodeCountX, grids[1].nodeCountY],counter)
		#plotParticles(problem.output,grid,counter)
		# se,ke    = computeEnergies(solids)
		# push!(problem.kinEnergy,ke)
		# push!(problem.strEnergy,se)
		# push!(problem.recordTime,t)
		#[compute(fix,t) for fix in problem.fixes]
	end
    t       += dtime
    counter += 1
  end # end of time loop
  #@printf("Solving done \n")
  #[closeFile(fix) for fix in problem.fixes]
end # end solve()

######################################################################
# Update Stress Last: FEM way
######################################################################
function solve_explicit_dynamics_2D_TL(grid,solids,basis,alg::FEM,output,Tf,dtime)
    t       = 0.
    counter = 0

    Identity = SMatrix{2,2}(1, 0, 0, 1)

	solidCount     = length(solids)
    interval       = output.interval

	linBasis       = LinearBasis()
	nearPointsLin  = [0,0,0,0]

    # pre_allocating arrays for temporary variables

 	Fdot          = zeros(Float64,2,2)
	D             = zeros(Float64,2,2)

	nearPoints,funcs, ders = initialise(basis,grids[1])

	if ( typeof(basis) <: CPDIQ4Basis )
		nearPointsLin    = [0, 0, 0, 0]
		funcsLin         = [0., 0., 0., 0.]
	end


    if ( typeof(output) <: PyPlotOutput )
	  pyFig_RealTime = PyPlot.figure(output.figTitle,
                               figsize=output.figSize, edgecolor="white", facecolor="white")
    end

    # compute once nodal mass
	# and probably shape functions ... and stored
	# and color field to recognise surface particles
	@inbounds for s = 1:solidCount
  	  solid  = solids[s]
	  grid   = grids[s]
	  # detect surface particles
	  mpoints = doCellParticleInteraction(solid,grid)
  	  # deformable solids only
  	  if solid.rigid continue end
  	  XX        = solid.pos0
  	  mm        = solid.mass
	  #rad       = solid.radius
	  vol       = solid.volumeInitial
  	  color     = solid.color
	  nodalMass = grid.mass
	  rp        = 1.0 / sqrt(grid.dx^2+grid.dy^2)
  	  @inbounds for ip = 1:solid.parCount
  		  support   = getShapeAndGradientTL(nearPoints,funcs,ders,ip, grid, solid, basis)
  		  fMass     = mm[ip]
		  #rad[ip]   = 0.5*cbrt(vol[ip])  # radius of the sphere of particle 'ip'
  		  #println(nearPoints)
  		  @inbounds for i = 1:support
  			  in              = nearPoints[i]; # index of node 'i'
  			  nodalMass[in]  += funcs[i] * fMass
  		  end
		  ################################################
		  # find surface particles by computing the sum kernel
		  xp = XX[ip][1]
		  yp = XX[ip][2]
		  e  = floor(Int64,(xp-grid.xmin)/grid.dx) + 1 +
		                   (grid.nodeCountX-1)*floor(Int64,(yp-grid.ymin)/grid.dy)
		  neighbors=getNeighbors(e, grid)
		  for ie=1:length(neighbors)
			  elemId = neighbors[ie]
			  pars   = mpoints[elemId]
			  #println(pars)
			  for iq = 1:length(pars)
				  q  = pars[iq]
				  xq = XX[q][1]
				  yq = XX[q][2]
				  dist = sqrt( (xp-xq)^2 +(yp-yq)^2 ) * rp
				  #println(dist)
                  if dist <= 1.
					  kernel     = 1-3*dist^2+2*dist^3
					  color[ip] += kernel
				  end
			  end
		  end
  	  end
    end

  while t < Tf

    @printf("Solving step: %d %f \n", counter, t)

    # ===========================================
    # reset grid data
    # ===========================================

	@inbounds for g=1:solidCount
		grid           = grids[g]
		#nodalMomentum0 = grid.momentum0
		nodalForce     = grid.force
	    @inbounds for i = 1:grid.nodeCount
		  #nodalMomentum0[i] = [0. 0.]
		  nodalForce[i]     = [0. 0.]
	    end
    end

    # ===========================================
    # particle to grid (deformable solids)
    # ===========================================

	@inbounds for s = 1:solidCount
		solid  = solids[s]
		grid   = grids[s]
		# deformable solids only
		if solid.rigid continue end

		nodalMomentum0 = grid.momentum0
		nodalForce     = grid.force

		XX     = solid.pos0
		mm     = solid.mass
		vv     = solid.velocity
		vol    = solid.volume
		stress = solid.stress
		F      = solid.deformationGradient
		vol0   = solid.volumeInitial

	  	@inbounds for ip = 1:solid.parCount
			support   = getShapeAndGradientTL(nearPoints,funcs,ders,ip, grid, solid, basis)
	        fVolume   = vol0[ip]
	        fMass     = mm[ip]
	        vp        = vv[ip]
	        sigma     = stress[ip]
			P         = det(F[ip])*sigma*inv(F[ip])'
			#body      = bodyforce(XX[ip],t)
			#println(nearPoints)
			@inbounds for i = 1:support
				in    = nearPoints[i]; # index of node 'i'
				Ni    = funcs[i]
				dNi   = ders[:,i]
				# momentum, internal force and external force and contact force
				#nodalMomentum0[in] += Ni * fMass * vp #+ vgrad*(grid.pos[in]-xx[ip]))
				nodalForce[in]     -=      fVolume * P     * dNi
				nodalForce[in]     +=      fMass   * body  *  Ni
			end
	  	end
	end

	# ===========================================
	# update grids (each grid for each solid)
	# ===========================================

	@inbounds for g=1:solidCount
		grid           = grids[g]
		nodalMomentum0 = grid.momentum0
		nodalMomentum  = grid.momentum
		nodalForce     = grid.force
		@inbounds for i=1:grid.nodeCount
			nodalMomentum[i] .= nodalMomentum0[i] + nodalForce[i] * dtime
	        # apply Dirichet boundary conditions
	        if grid.fixedXNodes[i] == 1
	        	nodalMomentum0[i][1]  = 0.
	        	nodalMomentum[i][1]   = 0.
	        end
	        if grid.fixedYNodes[i] == 1
	        	nodalMomentum0[i][2] = 0.
	        	nodalMomentum[i][2]  = 0.
	        end
		end
    end

    # ===========================================
    # grid to particle (deformable solids)
    # ===========================================

    @inbounds for s = 1:solidCount
	  	solid = solids[s]
		grid  = grids[s]

		nodalMass       = grid.mass
		nodalMomentum   = grid.momentum
		nodalMomentum0  = grid.momentum0

	  	xx    = solid.pos
	  	mm    = solid.mass
	  	vv    = solid.velocity
	  	vol   = solid.volume
	  	vol0  = solid.volumeInitial
	  	F     = solid.deformationGradient
	  	mat   = solid.mat
	  	stress = solid.stress
	  	strain = solid.strain
		vgrad  = solid.gradVelo

	  	@inbounds for ip = 1:solid.parCount
			support   = getShapeAndGradientTL(nearPoints,funcs,ders,ip, grid, solid, basis)
	        Fdot     .= 0.
			for i = 1:support
				in = nearPoints[i]; # index of node 'i'
				Ni = funcs[i]
				dNi= ders[:,i]
				mI = nodalMass[in]
			    if (mI > alg.tolerance)
					invM      = 1.0 / mI
					vI        = nodalMomentum[in] * invM
					vv[ip]   += Ni * (nodalMomentum[in] - nodalMomentum0[in]) * invM
					#vv[ip]   .= Ni * nodalMomentum[in] * invM, not compiled
					xx[ip]   += Ni * vI * dtime
					Fdot     += vI*dNi'
				end
			end
			F[ip]      += Fdot*dtime
			vel_grad    = Fdot*inv(F[ip])
			D           .= 0.5 * (vel_grad+vel_grad')
			strain[ip] += dtime * D
			J           = det(F[ip])
			vol[ip]     = J * vol0[ip]
			update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
	  	end
	end

    #### CPDI specific #################################
	### update corners
	if ( typeof(basis) <: CPDIQ4Basis )
		@inbounds for s = 1:solidCount
			solid = solids[s]
			if solid.rigid continue end
			@inbounds for c = 1:size(solid.nodes,1)
			  getShapeFuncs(nearPointsLin,funcsLin, solid.nodes[:,c], grid, solid, linBasis)
			  for i in 1:length(nearPointsLin)
				  in = nearPointsLin[i]; # index of node ‘i’
				  Ni = funcsLin[i]
				  #vI        = nodalMomentum2[in] / nodalMass[in]
				  if nodalMass[in] > alg.tolerance
					  solid.nodes[:,c]  += (Ni * nodalMomentum[in] / nodalMass[in]) * dtime
				  end
			  end
			end
		end
	end

	# ===========================================

	@inbounds for g=1:solidCount
		grid           = grids[g]
		nodalMomentum0 = grid.momentum0
		nodalMomentum  = grid.momentum
		@inbounds for i=1:grid.nodeCount
			nodalMomentum0[i] .= nodalMomentum[i]
		end
    end

	if (counter%interval == 0)
		plotParticles(output,solids,[grids[1].lx, grids[1].ly],
		              [grids[1].nodeCountX, grids[1].nodeCountY],counter)
		#plotParticles(problem.output,grid,counter)
		se,ke    = computeEnergies(solids)
		# push!(problem.kinEnergy,ke)
		# push!(problem.strEnergy,se)
		# push!(problem.recordTime,t)
		#[compute(fix,t) for fix in problem.fixes]
	end
    t       += dtime
    counter += 1
  end # end of time loop
  #@printf("Solving done \n")
  #[closeFile(fix) for fix in problem.fixes]
end # end solve()
