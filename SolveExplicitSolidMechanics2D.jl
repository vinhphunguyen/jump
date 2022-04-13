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
# USL and MUSL are provided, each function for each algorithm
#

using TimerOutputs
using Fix


######################################################################
# Modified Update Stress Last
######################################################################
function solve_explicit_dynamics_2D(grid,solids,basis,alg::MUSL,output,fixes,Tf,dtime)
	t             = 0.   # not t = 0 => type instability issue
    counter       = 0
    Identity      = UniformScaling(1.)
	solidCount    = length(solids)

	nodalMass      = grid.mass
	nodalMomentum  = grid.momentum
	nodalMomentum0 = grid.momentum0
	nodalMomentum2 = grid.momentum2
	nodalForce     = grid.force

	#vel_grad      = zeros(Float64,2,2)
	D             = SMatrix{2,2}(0., 0., 0., 0.) #zeros(Float64,2,2)
	vvp           = zeros(2)
	body          = zeros(2)

	alpha         = alg.alpha

	# allocate memory for grid basis and grads once
	show(basis)
	nearPoints,funcs, ders = initialise(grid,basis)

	if ( typeof(basis) <: CPDIQ4Basis )
		nearPoints = Vector{Int64}(undef,16)
		[nearPoints[i]=0 for i=1:16]
	    funcs = zeros(16)
		ders  = zeros(2,16)

		nearPointsLin    = [0, 0, 0, 0]
		funcsLin         = [0., 0., 0., 0.]
	end

    # if ( typeof(problem.output) <: PyPlotOutput )
	#   pyFig_RealTime = PyPlot.figure(problem.output.figTitle,
    #                            figsize=problem.output.figSize, edgecolor="white", facecolor="white")
    # end

	while t < Tf
	    #@printf(“Solving step: %f \n”, t)
	    # ===========================================
	    # reset grid data
	    # ===========================================
	     @inbounds for i = 1:grid.nodeCount
		  nodalMass[i]      = 0.
		  nodalMomentum0[i] = @SVector[0., 0.]
		  nodalMomentum2[i] = @SVector[0., 0.]
		  nodalForce[i]     = @SVector[0., 0.]
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
			F      = solid.deformationGradient
		  	@inbounds for ip = 1:solid.parCount
		        support   = getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid,basis)
		        fVolume   = vol[ip]
		        fMass     = mm[ip]
		        vp        = vv[ip]
		        sigma     = stress[ip]
		        Fm        = F[ip]
		        a         = Fm[1,1]
		        b         = Fm[1,2]
		        c         = Fm[2,1]
		        d         = Fm[2,2]
		        P         = sigma*SMatrix{2,2}(d, -c,-b, a) #zeros(Float64,2,2)[d -c;-b a] # convert to 1st PK stress
				#bodyforce(body,xx[ip],t)
					# println(nearPoints)
					# println(support)
				@inbounds for i = 1:support
					id    = nearPoints[i]; # index of node ‘i’
					Ni    = funcs[i]
					dNi   = @view ders[:,i]
					Nim   = Ni * fMass
					# mass, momentum, internal force and external force
					nodalMass[id]       +=  Nim
					nodalMomentum0[id]  +=  Nim * vp
					nodalForce[id]     -=      fVolume * P * dNi # (P[1,1] * dNi[1] + P[1,2] * dNi[2], P[2,1] * dNi[1] + P[2,2] * dNi[2])
	        #nodalForce[id] -= fVolume * @SVector[sigma[1,1] * dNi[1] + sigma[1,2] * dNi[2],sigma[2,1] * dNi[1] + sigma[2,2] * dNi[2]]
				end
		  	end
		end

		# ===========================================
		# update grid
		# ===========================================
		@inbounds for i=1:grid.nodeCount
		    nodalMomentum[i] = nodalMomentum0[i] + nodalForce[i] * dtime
	        # apply Dirichet boundary conditions
	       
        fixed_dirs       = @view grid.fixedNodes[:,i]
        if fixed_dirs[1] == 1
			    nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,1)
   		    nodalMomentum[i]  = setindex(nodalMomentum[i], 0.,1)
        end
        if fixed_dirs[2] == 1
			    nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,2)
			    nodalMomentum[i]  = setindex(nodalMomentum[i], 0.,2)
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
			  	@inbounds @simd for ip = 1:solid.parCount
					vp0     = vv[ip]
					mp      = mm[ip]
					vx     = 0.
					vy     = 0.
			        support = getShapeFunctions(nearPoints,funcs,ip, grid, solid,basis)

					for i =1:support
						in = nearPoints[i]; # index of node ‘i’
						Ni = funcs[i]
						mI = nodalMass[in]
						if ( mI > 0 )
							mii = Ni / mI
							vx += mii * (nodalMomentum[in][1] - alpha*nodalMomentum0[in][1])
							vy += mii * (nodalMomentum[in][2] - alpha*nodalMomentum0[in][2])
						end
					end
					vv[ip] = setindex(vv[ip],alpha*vp0[1] + vx,1)
					vv[ip] = setindex(vv[ip],alpha*vp0[2] + vy,2)
					# mapping the updated particle vel back to the node
					for i in 1:support
						in = nearPoints[i] # index of node ‘i’
					    nodalMomentum2[in]  += funcs[i] * mp * vv[ip]
					end
			  	end
	    end
	    # # apply Dirichet boundary conditions
	    @inbounds @simd for i = 1:grid.nodeCount
	       # apply Dirichet boundary conditions
        fixed_dirs       = @view grid.fixedNodes[:,i]
        if fixed_dirs[1] == 1
   		    nodalMomentum2[i]  = setindex(nodalMomentum2[i], 0.,1)
        end
        if fixed_dirs[2] == 1
			nodalMomentum2[i]  = setindex(nodalMomentum2[i], 0.,2)
        end

	     
	    end

	    # ===========================================
	    # particle to grid 2:
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
		  	@inbounds @simd for ip = 1:solid.parCount
				support   = getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid,basis)
		        #vel_grad .= 0. #zeros(Float64,2,2)
				vel_grad = SMatrix{2,2}(0., 0., 0., 0.)
				xxp = xx[ip]
				for i = 1:support
					in = nearPoints[i]; # index of node ‘i’
					Ni = funcs[i]
					dNi= @view ders[:,i]
					m  = nodalMass[in]
					if ( m > 0.)
						vI         = nodalMomentum2[in] /m
					    #xxp       += (Ni * nodalMomentum[in]/m) * dtime	
					    xxp       += (Ni * vI) * dtime	
				        vel_grad  += SMatrix{2,2}(dNi[1]*vI[1], dNi[1]*vI[2],
   										          dNi[2]*vI[1], dNi[2]*vI[2])
				    end
				end
				xx[ip]      = xxp
	            D           = 0.5 * (vel_grad + vel_grad')
	            strain[ip]  += dtime * D
				F[ip]       *= (Identity + vel_grad*dtime)
				J            = det(F[ip])
				if ( J < 0. )
					@printf("Troubled particle: %f %f \n", xx[ip][1], xx[ip][2])
					println(F[ip])
					@error("J is negative\n")
				end
				vol[ip]     = J * vol0[ip]
				#@timeit "3"  update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
				update_stress!(stress[ip],mat,strain[ip],D,F[ip],J,ip,dtime)
		  	end
		end

		# if CPDI, do extra thing here

		if ( typeof(basis) <: CPDIQ4Basis )
			@inbounds for s = 1:solidCount
				corner_coords = solids[s].nodes
				@inbounds for c = 1:length(corner_coords)
				  #xc = corner_coords[c]
				  #println(xc)
				  getShapeFuncs(nearPointsLin,funcsLin, corner_coords[c], grid, solids[s], LinearBasis())
				  for i = 1:length(nearPointsLin)
					  in = nearPointsLin[i]; # index of node ‘i’
					  Ni = funcsLin[i]
					  mI = nodalMass[in]
					  if mI > 0.
				        #corner_coords[c] += (Ni * nodalMomentum[in] / mI) * dtime
				        corner_coords[c] += (Ni * nodalMomentum2[in] / mI) * dtime
					  end
				  end
				  # println(xc)
				  # println(corner_coords[c])
				  # println("haha\n\n")
			    end
		    end
	    end

		if (counter%output.interval == 0)
			plotParticles_2D(output,solids,[grid.lx, grid.ly],
			             [grid.nodeCountX, grid.nodeCountY],counter)
			compute(fixes,t)
	    end

        t       += dtime
        counter += 1
    end
	closeFile(fixes)
end

######################################################################
# Update Stress Last
######################################################################
function solve_explicit_dynamics_2D(grid,solids,basis,alg::USL,output,fixes,Tf,dtime)
    t       = 0.
    counter = 0

    Identity       = UniformScaling(1.)
	solidCount     = length(solids)
	nodalMass      = grid.mass
	nodalMomentum0 = grid.momentum0
	nodalMomentum  = grid.momentum
	nodalForce     = grid.force

    D              = SMatrix{2,2}(0., 0., 0., 0.) #zeros(Float64,2,2)
	linBasis       = LinearBasis()
	nearPointsLin  = [0,0,0,0]

    # pre_allocating arrays for temporary variable

	nearPoints,funcs, ders = initialise(grid,basis)

	if ( typeof(basis) <: CPDIQ4Basis )
		nearPointsLin    = [0, 0, 0, 0]
		funcsLin         = [0., 0., 0., 0.]
	end


  while t < Tf

    #@printf("Solving step: %d%f \n", counter, t)

    # ===========================================
    # reset grid data
    # ===========================================

    @inbounds for i = 1:grid.nodeCount
	  nodalMass[i]      = 0.
	  nodalMomentum0[i] =  @SVector [0., 0.]
	  nodalForce[i]     =  @SVector [0., 0.]
    end

    # ===========================================
    # particle to grid (deformable solids)
    # ===========================================

	@inbounds for s = 1:solidCount
		solid  = solids[s]
		# deformable solids only
		if solid.rigid continue end
		xx     = solid.pos
		mm     = solid.mass
		vv     = solid.velocity
		vol    = solid.volume
		stress = solid.stress

	  	@inbounds for ip = 1:solid.parCount
	        #getShapeAndGradient(nearPoints,funcs,ders,xx[ip], grid)
			support   = getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid,basis)
	        fVolume   = vol[ip]
	        fMass     = mm[ip]
	        vp        = vv[ip]
	        sigma     = stress[ip]
			#body      = problem.bodyforce(xx[ip],t)
			#println(nearPoints)
			@inbounds for i = 1:support
				in    = nearPoints[i]; # index of node 'i'
				Ni    = funcs[i]
				dNi   = @view ders[:,i]
				Nim   = Ni * fMass
				# mass, momentum, internal force and external force
				nodalMass[in]      += Nim
				nodalMomentum0[in] += Nim * vp #+ vgrad*(grid.pos[in]-xx[ip]))
				nodalForce[in]     -= fVolume * @SVector[sigma[1,1] * dNi[1] + sigma[1,2] * dNi[2],
													     sigma[2,1] * dNi[1] + sigma[2,2] * dNi[2]]
				#nodalForce[in]     +=      fMass   * body  *  Ni
			end
	  	end
	end

	# ===========================================
	# update grid
	# ===========================================

	@inbounds for i=1:grid.nodeCount
		nodalMomentum[i] = nodalMomentum0[i] + nodalForce[i] * dtime
        # apply Dirichet boundary conditions
          # apply Dirichet boundary conditions
        fixed_dirs       = @view grid.fixedNodes[:,i]
        if fixed_dirs[1] == 1
			nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,1)
   		    nodalMomentum[i]  = setindex(nodalMomentum[i], 0.,1)
        end
        if fixed_dirs[2] == 1
			nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,2)
			nodalMomentum[i]  = setindex(nodalMomentum[i], 0.,2)
        end
	end

	# ===========================================
	# particle to grid (rigid solids)
	# ===========================================

	@inbounds for s = 1:solidCount
		solid  = solids[s]
		# deformable solids only
		if !solid.rigid continue end
		xx     = solid.pos
		vex     = solid.mat.vx
		vey     = solid.mat.vy
		@inbounds for ip = 1:solid.parCount
			getAdjacentGridPoints(nearPointsLin,xx[ip],grid,linBasis)
#			println(nearPoints)
			@inbounds for i = 1:4
				id                  = nearPointsLin[i]; # index of node 'i'
				mi                  = nodalMass[id]
				if solid.mat.fixed
				  nodalMomentum[id]   = setindex(nodalMomentum[id],0.,1)
				  nodalMomentum[id]   = setindex(nodalMomentum[id],0.,2)
  				  nodalMomentum0[id]  = setindex(nodalMomentum0[id],0.,1)
  				  nodalMomentum0[id]  = setindex(nodalMomentum0[id],0.,2)
			    else
					if vex != 0.
					  nodalMomentum[id]   = setindex(nodalMomentum[id], mi*vex,1)
					  nodalMomentum0[id]  = setindex(nodalMomentum0[id],mi*vex,1)
				    end
					if vey != 0.
					  nodalMomentum[id]   = setindex(nodalMomentum[id], mi*vey,2)
					  nodalMomentum0[id]  = setindex(nodalMomentum0[id],mi*vey,2)
					end
				end
				#println(nodalMomentum)
			end
		end
	end

    # ===========================================
    # grid to particle (deformable solids)
    # ===========================================

    @inbounds for s = 1:solidCount
		# only deformable solids here
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

	  	@inbounds for ip = 1:solid.parCount
			support   = getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid,basis)
	        vel_grad  = SMatrix{2,2}(0., 0., 0., 0.)
			vvp       = vv[ip]
			xxp       = xx[ip]
			for i = 1:support
				in = nearPoints[i]; # index of node 'i'
				Ni = funcs[i]
				dNi= @view ders[:,i]
				mI = nodalMass[in]
			    if (mI > alg.tolerance)
					invM       = 1.0 / mI
					vI         = nodalMomentum[in] * invM
					vvp       += Ni * (nodalMomentum[in] - nodalMomentum0[in]) * invM
					xxp       += Ni * vI * dtime
			        vel_grad  += SMatrix{2,2}(dNi[1]*vI[1], dNi[1]*vI[2],
   										      dNi[2]*vI[1], dNi[2]*vI[2])
		   		end
		   	 end
			 vv[ip]      = vvp
			 xx[ip]      = xxp
		   	 D           = 0.5 * (vel_grad + vel_grad')
		   	 strain[ip]  += dtime * D
		   	 F[ip]       *= (Identity + vel_grad*dtime)
		   	 J            = det(F[ip])
		   	 # if ( J < 0. )
		   		#  @printf("Troubled particle: %f %f \n", xx[ip][1], xx[ip][2])
		   		#  println(F[ip])
		   		#  @error("J is negative\n")
		   	 # end
	   	 vol[ip]     = J * vol0[ip]
	   	 #@timeit "3" update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
	   	 update_stress!(stress[ip],mat,strain[ip],D,F[ip],J,ip,dtime)
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
	  	vx    = solid.mat.vx
	  	vy    = solid.mat.vy
		#println(ve)
        @inbounds for ip = 1:solid.parCount
	      xx[ip]   += dtime * @SVector [vx,vy]
	    end
	end

    #### CPDI specific #################################
	### update corners

	if ( typeof(basis) <: CPDIQ4Basis )
		@inbounds for s = 1:solidCount
			corner_coords = solids[s].nodes
			@inbounds for c = 1:length(corner_coords)
			  xc = corner_coords[c]
			  getShapeFuncs(nearPointsLin,funcsLin, xc, grid, solids[s], linBasis)
			  for i = 1:length(nearPointsLin)
				  in = nearPointsLin[i]; # index of node ‘i’
				  Ni = funcsLin[i]
				  mI = nodalMass[in]
				  if mI > 0.
					xc   += (Ni * nodalMomentum[in] / mI) * dtime
				  end
			  end
			end
		end
	end

	if (counter%output.interval == 0)
		plotParticles_2D(output,solids,[grid.lx, grid.ly],
					 [grid.nodeCountX, grid.nodeCountY],counter)
		compute(fixes,t)
	end

    t       += dtime
    counter += 1
  end # end of time loop
  #@printf("Solving done \n")
  closeFile(fixes)
end # end solve()
