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



######################################################################
# Update Stress Last
######################################################################
function solve_explicit_dynamics_femp_2D(grid,solids,basis,bodyforce,alg::USL,output,fixes,Tf,dtime)
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

    dNdx      = zeros(2,4)
    vel_grad  = SMatrix{2,2}(0., 0., 0., 0.)
    N         = zeros(4)#@SVector [0,0,0,0]

    # compute nodal mass (only once)

     gpCoords  = zeros(2,1)
    weights   = [4.] #ones(1)
    # gpCoords[1,1] = -0.5773502691896257; gpCoords[2,1] = -0.5773502691896257;
    # gpCoords[1,2] =  0.5773502691896257; gpCoords[2,2] = -0.5773502691896257;
    # gpCoords[1,3] =  0.5773502691896257; gpCoords[2,3] =  0.5773502691896257;
    # gpCoords[1,4] = -0.5773502691896257; gpCoords[2,4] =  0.5773502691896257;

    noGP = 1
    gpCoords=[0.,0.]

    @inbounds for s = 1:solidCount
		solid  = solids[s]
		xx     = solid.pos
		elems  = solid.elems

	  	@inbounds for ip = 1:solid.parCount
			elemNodes  =  @view elems[ip,:]  
			elemNodes0 =        elems[ip,:]  
			coords     =  @view xx[elemNodes]
	    
			@inbounds for gp = 1:noGP
			  xieta = @view gpCoords[:,gp]
			  detJ  = lagrange_basis!(N, Quad4(), xieta, coords)
			  if detJ < 0.
			  	elemNodes[2] = elemNodes0[4]			  	
			  	elemNodes[4] = elemNodes0[2]
			  	#println(N)
			  end
			end
	  	end
    end

    @inbounds for s = 1:solidCount
		solid  = solids[s]
		xx     = solid.pos
		mm     = solid.mass  # to be updated here
		elems  = solid.elems
		mat    = solid.mat
		vol    = solid.volume

	  	@inbounds for ip = 1:solid.parCount
			elemNodes  =  @view elems[ip,:]  
			elemNodes0 =        elems[ip,:]  
			coords     =  @view xx[elemNodes]
	    
			@inbounds for gp = 1:noGP
			  xieta = @view gpCoords[:,gp]
			  detJ  = lagrange_basis!(N, Quad4(), xieta, coords)
			  if detJ < 0.
			  	elemNodes[2] = elemNodes0[4]			  	
			  	elemNodes[4] = elemNodes0[2]
			  	println(N)
			  end
			  	
			  #println(elemNodes)	
			  for i=1:length(elemNodes)
			  	id      = elemNodes[i]
			  	mm[id]  += mat.density*N[i]*detJ*weights[gp]
			  	vol[ip] += detJ*weights[gp]
			  end
			end
	  	end
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
	
		xx     = solid.pos
		mm     = solid.mass
		vv     = solid.velocity		
		fint   = solid.fint
		du     = solid.dU

	  	@inbounds for ip = 1:solid.nodeCount
	        #getShapeAndGradient(nearPoints,funcs,ders,xx[ip], grid)
			support   = getShapeFunctions(nearPoints,funcs,ip, grid, solid,basis)
	        
	        fMass     = mm[ip]
	        vp        = vv[ip]
	        fp        = fint[ip]
			#body      = problem.bodyforce(xx[ip],t)
			#println(nearPoints)
			@inbounds for i = 1:support
				in    = nearPoints[i]; # index of node 'i'
				Ni    = funcs[i]				
				Nim   = Ni * fMass
				# mass, momentum, internal force and external force
				nodalMass[in]      += Nim
				nodalMomentum0[in] += Nim * vp 
				nodalForce[in]     -= Ni  * fp
				#nodalForce[in]     +=      fMass   * body  *  Ni
			end
			#du[ip] = @SVector [0., 0.]
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
	end

	
    # ====================================================================
    # grid to particle (deformable solids): update particle velocity/pos
    # ====================================================================

    @inbounds for s = 1:solidCount
		# only deformable solids here
	  	solid = solids[s]

	  	xx    = solid.pos
	  	mm    = solid.mass
	  	vv    = solid.velocity
	  	du    = solid.dU
	  	
	  	@inbounds for ip = 1:solid.nodeCount
			support   = getShapeFunctions(nearPoints,funcs,ip, grid, solid, basis)	        
			vvp       = vv[ip]
			xxp       = xx[ip]
			dup       = du[ip]
			for i = 1:support
				in = nearPoints[i]; # index of node 'i'
				Ni = funcs[i]
				mI = nodalMass[in]
			    if (mI > alg.tolerance)
					invM       = 1.0 / mI
					vI         = nodalMomentum[in] * invM
					vvp       += Ni * (nodalMomentum[in] - nodalMomentum0[in]) * invM
					xxp       += Ni * vI * dtime				
					dup       += Ni * vI * dtime				
		   		end
		   	end
			vv[ip]      = vvp
			xx[ip]      = xxp	
			du[ip]      = dup
	    end
	end


    # ====================================================================
    #  update particle internal forces
    # ====================================================================

    @inbounds for s = 1:solidCount
		# only deformable solids here
	  	solid = solids[s]

	  	xx     = solid.pos
	  	mm     = solid.mass
	  	du     = solid.dU
	  	F      = solid.deformationGradient
	  	mat    = solid.mat
	  	stress = solid.stress
	  	strain = solid.strain
	  	elems  = solid.elems
	  	fint   = solid.fint
	  	vol    = solid.volume
	  	for ip=1:solid.nodeCount 
	  		fint[ip] = @SVector [0., 0.]
	  	end
        # loop over solid elements, not solid nodes
	  	@inbounds for ip = 1:solid.parCount
			elemNodes =  @view elems[ip,:]  
			coords    =  @view xx[elemNodes]
	        vel_grad  =  SMatrix{2,2}(0., 0., 0., 0.)
			
			detJ      = lagrange_basis_derivatives!(dNdx, Quad4(), [0. 0.], coords)
			vol[ip]   = detJ * 4 # weight of GP = 4
			#println(dNdx)
			#println(sum(dNdx, dims=2))
			for i = 1:length(elemNodes)
				in         = elemNodes[i]; # index of node 'i'
			    dNi        = @view dNdx[:,i]			
				vI         = du[in]
				vel_grad  += SMatrix{2,2}(dNi[1]*vI[1], dNi[2]*vI[1],
   										  dNi[1]*vI[2], dNi[2]*vI[2])
		   	end
	
		   	D           = 0.5 * (vel_grad + vel_grad')
		   	strain[ip]  = D
		   	F[ip]        = inv(Identity - vel_grad)
		   	J            = det(F[ip])
		   	#println(strain[ip])
	   	     #@timeit "3" update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
	   	    update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)

            sigma = stress[ip]
            # compute nodal internal force fint
		    for i = 1:length(elemNodes)
				in  = elemNodes[i]; # index of node 'i'
			    dNi = @view dNdx[:,i]			
	   	        fint[in]     += detJ * 4 * @SVector[sigma[1,1] * dNi[1] + sigma[1,2] * dNi[2],
											        sigma[2,1] * dNi[1] + sigma[2,2] * dNi[2]]
            end
	   end
	end

	if (counter%output.interval == 0)
		plotParticles_2D(output,solids,counter)
		compute_femp(fixes,t)
	end

    t       += dtime
    counter += 1
  end # end of time loop
  #@printf("Solving done \n")
  closeFile(fixes)
end # end solve()


######################################################################
# Update Stress Last, TLFEM for internal force
######################################################################
function solve_explicit_dynamics_femp_2D(grid,solids,basis,body,alg::TLFEM,output,fixes,Tf,dtime)
    t       = 0.
    counter = 0

    Identity       = UniformScaling(1.)
	solidCount     = length(solids)
	nodalMass      = grid.mass
	nodalMomentum0 = grid.momentum0
	nodalMomentum  = grid.momentum
	nodalForce     = grid.force


    # pre_allocating arrays for temporary variable
   
	nearPoints,funcs, ders = initialise(grid,basis)

    dNdx      = zeros(2,4)
    vel_grad  = SMatrix{2,2}(0., 0., 0., 0.)
    D         = SMatrix{2,2}(0., 0., 0., 0.) #zeros(Float64,2,2)
    g         = [0.,0.]
    N         = zeros(4)#@SVector [0,0,0,0]

    # compute nodal mass (only once)

    gpCoords  = zeros(2,1)
    weights   = [4.] #ones(1)
    # gpCoords[1,1] = -0.5773502691896257; gpCoords[2,1] = -0.5773502691896257;
    # gpCoords[1,2] =  0.5773502691896257; gpCoords[2,2] = -0.5773502691896257;
    # gpCoords[1,3] =  0.5773502691896257; gpCoords[2,3] =  0.5773502691896257;
    # gpCoords[1,4] = -0.5773502691896257; gpCoords[2,4] =  0.5773502691896257;

    noGP = 1
    gpCoords=[0.,0.]

    @inbounds for s = 1:solidCount
		solid  = solids[s]
		xx     = solid.pos
		elems  = solid.elems

	  	@inbounds for ip = 1:solid.parCount
			elemNodes  =  @view elems[ip,:]  
			elemNodes0 =        elems[ip,:]  
			coords     =  @view xx[elemNodes]
	    
			@inbounds for gp = 1:noGP
			  xieta = @view gpCoords[:,gp]
			  detJ  = lagrange_basis!(N, Quad4(), xieta, coords)
			  if detJ < 0.
			  	elemNodes[2] = elemNodes0[4]			  	
			  	elemNodes[4] = elemNodes0[2]
			  	#println(N)
			  end
			end
	  	end
    end

    @inbounds for s = 1:solidCount
		solid  = solids[s]
		xx     = solid.pos
		mm     = solid.mass  # to be updated here
		elems  = solid.elems
		mat    = solid.mat
		vol    = solid.volume

	  	@inbounds for ip = 1:solid.parCount
			elemNodes  =  @view elems[ip,:]  
			elemNodes0 =        elems[ip,:]  
			coords     =  @view xx[elemNodes]
	    
			@inbounds for gp = 1:noGP
			  xieta = @view gpCoords[:,gp]
			  detJ  = lagrange_basis!(N, Quad4(), xieta, coords)
			  if detJ < 0.
			  	elemNodes[2] = elemNodes0[4]			  	
			  	elemNodes[4] = elemNodes0[2]
			  	println(N)
			  end
			  	
			  #println(elemNodes)	
			  for i=1:length(elemNodes)
			  	id      = elemNodes[i]
			  	mm[id]  += mat.density*N[i]*detJ*weights[gp]
			  	vol[ip] += detJ*weights[gp]
			  end
			end
	  	end
    end

  while t < Tf

    @printf("Solving step: %d %f \n", counter, t)

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
	
		xx     = solid.pos
		mm     = solid.mass
		vv     = solid.velocity		
		fint   = solid.fint
		fbody  = solid.fbody
		du     = solid.dU

	  	@inbounds for ip = 1:solid.nodeCount
	        #getShapeAndGradient(nearPoints,funcs,ders,xx[ip], grid)
			support   = getShapeFunctions(nearPoints,funcs,ip, grid, solid,basis)
	        
	        fMass     = mm[ip]
	        vp        = vv[ip]
	        fp        = fint[ip]
	        fb        = fbody[ip]
			
			#println(funcs)
			@inbounds for i = 1:support
				in    = nearPoints[i]; # index of node 'i'
				Ni    = funcs[i]				
				Nim   = Ni * fMass
				# mass, momentum, internal force and external force
				nodalMass[in]      += Nim
				nodalMomentum0[in] += Nim * vp 
				nodalForce[in]     -= Ni  * fp
				nodalForce[in]     += Ni  * fb
			end
			#du[ip] = @SVector [0., 0.]
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
	end

	
    # ====================================================================
    # grid to particle (deformable solids): update particle velocity/pos
    # ====================================================================

    @inbounds for s = 1:solidCount
		# only deformable solids here
	  	solid = solids[s]

	  	xx    = solid.pos
	  	mm    = solid.mass
	  	vv    = solid.velocity
	  	du    = solid.dU
	    fixX  = solid.fixedNodesX
	  	fixY  = solid.fixedNodesY
	  	
	  	
	  	@inbounds for ip = 1:solid.nodeCount
			support   = getShapeFunctions(nearPoints,funcs,ip, grid, solid, basis)	        
			vvp       = vv[ip]
			xxp       = xx[ip]
			dup       = du[ip]
			#vvt = [0.,0.]
			for i = 1:support
				in = nearPoints[i]; # index of node 'i'
				Ni = funcs[i]
				mI = nodalMass[in]
			    if (mI > alg.tolerance)
					invM       = 1.0 / mI
					vI         = nodalMomentum[in] * invM
					#vvt        += Ni * vI  => too much dissipation
					vvp       += Ni * (nodalMomentum[in] - nodalMomentum0[in]) * invM
					xxp       += Ni * vI * dtime				
					dup       += Ni * vI * dtime				
		   		end
		   	end
			vv[ip]      = vvp
			xx[ip]      = xxp	
			du[ip]      = dup
		    # Dirichlet BCs on the mesh
			if fixX[ip] == 1 
				vv[ip] = setindex(vv[ip],0.,1)
				du[ip] = setindex(du[ip],0.,1)
				xx[ip] = setindex(xx[ip],solid.pos0[ip][1],1)
			end
				# Dirichlet BCs on the mesh
			if fixY[ip] == 1 
				vv[ip] = setindex(vv[ip],0.,2)
				du[ip] = setindex(du[ip],0.,2)
				xx[ip] = setindex(xx[ip],solid.pos0[ip][2],2)
			end
	    end
	end


    # ====================================================================
    #  update particle internal forces
    # ====================================================================
t       += dtime
    @inbounds for s = 1:solidCount
		# only deformable solids here
	  	solid = solids[s]

	  	XX     = solid.pos0
	  	mm     = solid.mass
	  	du     = solid.dU
	  	velo    = solid.velocity
	  	F      = solid.deformationGradient
	  	mat    = solid.mat
	  	stress = solid.stress
	  	strain = solid.strain
	  	elems  = solid.elems
	  	fint   = solid.fint
	  	fbody  = solid.fbody
	  	vol    = solid.volume
	  	for ip=1:solid.nodeCount 
	  		fint[ip]  = @SVector [0., 0.]
	  		fbody[ip] = @SVector [0., 0.]
	  	end
        # loop over solid elements, not solid nodes
	  	@inbounds for ip = 1:solid.parCount
			elemNodes =  @view elems[ip,:]  
			coords    =  @view XX[elemNodes]
	        vel_grad  =  SMatrix{2,2}(0., 0., 0., 0.)
			
			detJ      = lagrange_basis_derivatives!(N, dNdx, Quad4(), [0. 0.], coords)
			w         = detJ * 4
			vol[ip]   = w
			#println(dNdx)
			#println(sum(dNdx, dims=2))
			for i = 1:length(elemNodes)
				in         = elemNodes[i]; # index of node 'i'
			    dNi        = @view dNdx[:,i]			
				vI         = du[in]
				#vI         = velo[in]
				vel_grad  += SMatrix{2,2}(dNi[1]*vI[1], dNi[1]*vI[2],
   										  dNi[2]*vI[1], dNi[2]*vI[2])
		   	end
			   	
		   	F[ip]        = Identity + vel_grad
		   	#F[ip]        += dtime * vel_grad
		   	# D            = 0.5 * (vel_grad + vel_grad')
		   	strain[ip]   = 0.5 * (vel_grad + vel_grad' + vel_grad * vel_grad')

		   	# L              = vel_grad * inv(F[ip])
		   	# D              = 0.5 * (L + L')
		   	# strain[ip]    += dtime * D

		   	J            = det(F[ip])
		   	if ( J < 0. )
					@printf("Troubled particle: %d %f %f \n", s, XX[ip][1], XX[ip][2])
					#println(F[ip])
					#@error("J is negative\n")
			end
		   	#println(strain[ip])
	   	     #@timeit "3" update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
	   	    update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)

            xcenter = @SVector[0.25*(coords[1][1]+coords[2][1]+coords[3][1]+coords[4][1]), 
                               0.25*(coords[1][2]+coords[2][2]+coords[3][2]+coords[4][2])]
            body(g,xcenter,t)  
            sigma = stress[ip]
            P     = J*sigma*inv(F[ip])'  # convert to 1st Piola Kirchoof stress
            # compute nodal internal force fint
		    for i = 1:length(elemNodes)
				in  = elemNodes[i]; # index of node 'i'
			    dNi = @view dNdx[:,i]			
	   	        fint[in]     += w * @SVector[P[1,1] * dNi[1] + P[1,2] * dNi[2],
											 P[2,1] * dNi[1] + P[2,2] * dNi[2]]
                fbody[in]    += w*mat.density*0.25*g															        
                #fbody[in]    += w*mat.density*N[i]*g															        
            end
	   end
	end

	if (counter%output.interval == 0)
		plotParticles_2D(output,solids,counter)
		compute_femp(fixes,t)
	end

#    t       += dtime
    counter += 1
  end # end of time loop
  #@printf("Solving done \n")
  closeFile(fixes)
end # end solve()



######################################################################
# Update Stress Last, TLFEM for internal force, full Gauss quadrature
######################################################################
function solve_explicit_dynamics_femp_2D(grid,solids,basis,body,alg::TLFEMFull,output,fixes,Tf,dtime)
    t       = 0.
    counter = 0

    Identity       = UniformScaling(1.)
	solidCount     = length(solids)
	nodalMass      = grid.mass
	nodalMomentum0 = grid.momentum0
	nodalMomentum  = grid.momentum
	nodalForce     = grid.force


    # pre_allocating arrays for temporary variable
   
	nearPoints,funcs, ders = initialise(grid,basis)

    dNdx      = zeros(2,4)
    vel_grad  = SMatrix{2,2}(0., 0., 0., 0.)
    D         = SMatrix{2,2}(0., 0., 0., 0.) #zeros(Float64,2,2)
    g         = [0.,0.]
    N         = zeros(4)#@SVector [0,0,0,0]

    # compute nodal mass (only once)

    gpCoords  = zeros(2,4)
    weights   = ones(1)
    gpCoords[1,1] = -0.5773502691896257; gpCoords[2,1] = -0.5773502691896257;
    gpCoords[1,2] =  0.5773502691896257; gpCoords[2,2] = -0.5773502691896257;
    gpCoords[1,3] =  0.5773502691896257; gpCoords[2,3] =  0.5773502691896257;
    gpCoords[1,4] = -0.5773502691896257; gpCoords[2,4] =  0.5773502691896257;

    noGP = 4
    # gpCoords=[0.,0.]

    @inbounds for s = 1:solidCount
		solid  = solids[s]
		xx     = solid.pos
		elems  = solid.elems

	  	@inbounds for ip = 1:solid.parCount
			elemNodes  =  @view elems[ip,:]  
			elemNodes0 =        elems[ip,:]  
			coords     =  @view xx[elemNodes]
	    
			@inbounds for gp = 1:noGP
			  xieta = @view gpCoords[:,gp]
			  detJ  = lagrange_basis!(N, Quad4(), xieta, coords)
			  if detJ < 0.
			  	elemNodes[2] = elemNodes0[4]			  	
			  	elemNodes[4] = elemNodes0[2]
			  	#println(N)
			  end
			end
	  	end
    end

    @inbounds for s = 1:solidCount
		solid  = solids[s]
		xx     = solid.pos
		mm     = solid.mass  # to be updated here
		elems  = solid.elems
		mat    = solid.mat
		vol    = solid.volume

	  	@inbounds for ip = 1:solid.parCount
			elemNodes  =  @view elems[ip,:]  
			elemNodes0 =        elems[ip,:]  
			coords     =  @view xx[elemNodes]
	    
			@inbounds for gp = 1:noGP
			  xieta = @view gpCoords[:,gp]
			  detJ  = lagrange_basis!(N, Quad4(), xieta, coords)
			  if detJ < 0.
			  	elemNodes[2] = elemNodes0[4]			  	
			  	elemNodes[4] = elemNodes0[2]
			  	println(N)
			  end
			  	
			  #println(elemNodes)	
			  for i=1:length(elemNodes)
			  	id      = elemNodes[i]
			  	mm[id]  += mat.density*N[i]*detJ*weights[gp]
			  	vol[ip] += detJ*weights[gp]
			  end
			end
	  	end
    end

  stress1 = zeros(2,2)
  F1      = zeros(2,2)


  while t < Tf

    @printf("Solving step: %d %f \n", counter, t)

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
	
		xx     = solid.pos
		mm     = solid.mass
		vv     = solid.velocity		
		fint   = solid.fint
		fbody  = solid.fbody
		du     = solid.dU

	  	@inbounds for ip = 1:solid.nodeCount
	        #getShapeAndGradient(nearPoints,funcs,ders,xx[ip], grid)
			support   = getShapeFunctions(nearPoints,funcs,ip, grid, solid,basis)
	        
	        fMass     = mm[ip]
	        vp        = vv[ip]
	        fp        = fint[ip]
	        fb        = fbody[ip]
			#body      = problem.bodyforce(xx[ip],t)
			#println(nearPoints)
			@inbounds for i = 1:support
				in    = nearPoints[i]; # index of node 'i'
				Ni    = funcs[i]				
				Nim   = Ni * fMass
				# mass, momentum, internal force and external force
				nodalMass[in]      += Nim
				nodalMomentum0[in] += Nim * vp 
				nodalForce[in]     -= Ni  * fp
				nodalForce[in]     += Ni  * fb
			end
			#du[ip] = @SVector [0., 0.]
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
	end

	
    # ====================================================================
    # grid to particle (deformable solids): update particle velocity/pos
    # ====================================================================

    @inbounds for s = 1:solidCount
		# only deformable solids here
	  	solid = solids[s]

	  	xx    = solid.pos
	  	mm    = solid.mass
	  	vv    = solid.velocity
	  	du    = solid.dU
	    fixX  = solid.fixedNodesX
	  	fixY  = solid.fixedNodesY
	  	
	  	
	  	@inbounds for ip = 1:solid.nodeCount
			support   = getShapeFunctions(nearPoints,funcs,ip, grid, solid, basis)	        
			vvp       = vv[ip]
			xxp       = xx[ip]
			dup       = du[ip]
			for i = 1:support
				in = nearPoints[i]; # index of node 'i'
				Ni = funcs[i]
				mI = nodalMass[in]
			    if (mI > alg.tolerance)
					invM       = 1.0 / mI
					vI         = nodalMomentum[in] * invM
					vvp       += Ni * (nodalMomentum[in] - nodalMomentum0[in]) * invM
					xxp       += Ni * vI * dtime				
					dup       += Ni * vI * dtime				
		   		end
		   	end
			vv[ip]      = vvp
			xx[ip]      = xxp	
			du[ip]      = dup
		    # Dirichlet BCs on the mesh
			if fixX[ip] == 1 
				vv[ip] = setindex(vv[ip],0.,1)
				du[ip] = setindex(du[ip],0.,1)
				xx[ip] = setindex(xx[ip],solid.pos0[ip][1],1)
			end
				# Dirichlet BCs on the mesh
			if fixY[ip] == 1 
				vv[ip] = setindex(vv[ip],0.,2)
				du[ip] = setindex(du[ip],0.,2)
				xx[ip] = setindex(xx[ip],solid.pos0[ip][2],2)
			end
	    end
	end


    # ====================================================================
    #  update particle internal forces
    # ====================================================================
t       += dtime
    @inbounds for s = 1:solidCount
		# only deformable solids here
	  	solid = solids[s]

	  	XX     = solid.pos0
	  	mm     = solid.mass
	  	du     = solid.dU
	  	velo   = solid.velocity
	  	F      = solid.deformationGradient
	  	mat    = solid.mat
	  	stress = solid.stress
	  	strain = solid.strain
	  	elems  = solid.elems
	  	fint   = solid.fint
	  	fbody  = solid.fbody
	  	vol    = solid.volume
	  	for ip=1:solid.nodeCount 
	  		fint[ip]  = @SVector [0., 0.]
	  		fbody[ip] = @SVector [0., 0.]
	  	end
        # loop over solid elements, not solid nodes
	  	@inbounds for ip = 1:solid.parCount
			elemNodes =  @view elems[ip,:]  
			coords    =  @view XX[elemNodes]
	        vel_grad  =  SMatrix{2,2}(0., 0., 0., 0.)

	        @inbounds for gp = 1:noGP
			    xieta = @view gpCoords[:,gp]
			
			    detJ      = lagrange_basis_derivatives!(N, dNdx, Quad4(), xieta, coords)

			    w         = detJ * weights[gp]
		
				for i = 1:length(elemNodes)
					in         = elemNodes[i]; # index of node 'i'
				    dNi        = @view dNdx[:,i]			
					#vI         = du[in]
					vI         = velo[in]
					vel_grad  += SMatrix{2,2}(dNi[1]*vI[1], dNi[1]*vI[2],
	   										  dNi[2]*vI[1], dNi[2]*vI[2])
			   	end
		
			   	D            = 0.5 * (vel_grad + vel_grad')
			   	#strain[ip]   = D
			   	F1           = Identity + vel_grad
			   	J            = det(F1)
			   	if ( J < 0. )
						@printf("Troubled particle: %f %f \n", XX[ip][1], XX[ip][2])
						println(F1)
						@error("J is negative\n")
				end
			   	#println(strain[ip])
		   	     #@timeit "3" update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)

		   	    #update_stress!(stress,mat,strain[0],F,J,0)
# println(stress)
# println(F)
		   	    stress1    .= (1.0/J)*( mat.mu*(F1*F1'-UniformScaling(1.)) + mat.lambda*log(J)*UniformScaling(1.) )

	            xcenter = @SVector[(N[1]*coords[1][1]+N[2]*coords[2][1]+N[3]*coords[3][1]+N[4]*coords[4][1]), 
	                               (N[1]*coords[1][2]+N[2]*coords[2][2]+N[3]*coords[3][2]+N[4]*coords[4][2])]
	            body(g,xcenter,t)  
	            sigma = stress1
	            P     = J*sigma*inv(F1)'  # convert to 1st Piola Kirchoof stress
	            # compute nodal internal force fint
			    for i = 1:length(elemNodes)
					in  = elemNodes[i]; # index of node 'i'
				    dNi = @view dNdx[:,i]			
		   	        fint[in]     += w * @SVector[P[1,1] * dNi[1] + P[1,2] * dNi[2],
												        P[2,1] * dNi[1] + P[2,2] * dNi[2]]
	                fbody[in]    += w*mat.density*N[i]*g															        
	            end
	        end    
	   end
	end

	if (counter%output.interval == 0)
		plotParticles_2D(output,solids,counter)
		compute_femp(fixes,t)
	end

#    t       += dtime
    counter += 1
  end # end of time loop
  #@printf("Solving done \n")
  closeFile(fixes)
end # end solve()
