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
# Update Stress Last, TLFEM for internal force
######################################################################
function solve_explicit_dynamics_femp_2D_Contact(grid,solids,mats,basis,body,alg::TLFEM,output,fixes,data)
    Tf    = data["total_time"]
	dtime = data["dt"]         
	t     = data["time"]      
	fric  = data["friction"]

    counter = 0

    Identity       = UniformScaling(1.)
	solidCount     = length(solids)

    # this is body mass/momenta/forces/normals
    nodalMass      =  [copy(grid.mass)       for _ in 1:solidCount]
    nodalMomentum0 =  [copy(grid.momentum0)  for _ in 1:solidCount]
    nodalMomentum  =  [copy(grid.momentum)   for _ in 1:solidCount]
    nodalForce     =  [copy(grid.force)      for _ in 1:solidCount]
    nodalNormals   =  [copy(grid.force)      for _ in 1:solidCount]

    # this is system mass/momenta/force
    nodalMass_S      = grid.mass
	nodalMomentum0_S = grid.momentum0
	nodalMomentum_S  = grid.momentum
	nodalForce_S     = grid.force

	alpha            = alg.alpha

    nodeCount        = grid.nodeCount

    # for all grid nodes, store ids of solids that are in contact
    contact_solids   = Vector{Set{Int64}}(undef,nodeCount)



    # pre_allocating arrays for temporary variable
   
	nearPoints,funcs, ders = initialise(grid,basis)

    # for rigid bodies
	linBasis       = LinearBasis()
	nearPointsLin  = [0,0,0,0]

    
    vel_grad  = SMatrix{2,2}(0., 0., 0., 0.)
    D         = SMatrix{2,2}(0., 0., 0., 0.) #zeros(Float64,2,2)
    Fdot      = zeros(2,2)
    g         = [0.,0.]
    deltaVe1  = zeros(2)
    dVxNI     = zeros(2)
    omega     = zeros(2)
    nI        = zeros(2)
    nI_common = zeros(2)
   

    # compute nodal mass (only once)

    noGP        = 1
    nodePerElem = size(solids[1].elems,2)

	if nodePerElem == 3 		
		wgt       = .5
		weights   = [0.5]   # this is so dangerous!!! all books ay weight =1
		gpCoords  = [0.3333333333333,0.3333333333333]
		N         = zeros(3)#@SVector [0,0,0,0]
		dNdx      = zeros(2,3)
	end


	if nodePerElem == 4		
		wgt       = 4.
		weights   = [4.0]   # this is so dangerous!!! all books ay weight =1
		gpCoords  = [0.0,0.0]
		N         = zeros(4)#@SVector [0,0,0,0]
		dNdx      = zeros(2,4)
	end

    @inbounds for s = 1:solidCount
		solid  = solids[s]
		xx     = solid.pos
		elems  = solid.elems
		meshBasis = solid.basis

	  	@inbounds for ip = 1:solid.parCount
			elemNodes  =  @view elems[ip,:]  
			elemNodes0 =        elems[ip,:]  
			coords     =  @view xx[elemNodes]
	    
			@inbounds for gp = 1:noGP
			  xieta = @view gpCoords[:,gp]
			  detJ  = lagrange_basis!(N, meshBasis, xieta, coords)
			  if detJ < 0.
			  	if nodePerElem == 4
			  	    elemNodes[2] = elemNodes0[4]			  	
			  	    elemNodes[4] = elemNodes0[2]
			    elseif nodePerElem == 3
                    elemNodes[2] = elemNodes0[3]			  	
			  	    elemNodes[3] = elemNodes0[2]
			    end
			  	#println(N)
			  end
			end
	  	end
    end

    @inbounds for s = 1:solidCount
		solid  = solids[s]
		if typeof(mats[s]) <: RigidMaterial 
			compute_normals!(solid)
			continue 
		end
		xx     = solid.pos
		mm     = solid.mass  # to be updated here
		elems  = solid.elems
		rho    = mats[s].density
		vol    = solid.volume
		meshBasis = solid.basis

	  	@inbounds for ip = 1:solid.parCount
			elemNodes  =  @view elems[ip,:]  
			elemNodes0 =        elems[ip,:]  
			coords     =  @view xx[elemNodes]
	    
			@inbounds for gp = 1:noGP
			  xieta = @view gpCoords[:,gp]
			  detJ  = lagrange_basis!(N, meshBasis, xieta, coords)
			  if detJ < 0.
			  	println(N)
			  end
			  vol[ip] += detJ*weights[gp]
			  	
			  #println(elemNodes)	
			  for i=1:length(elemNodes)
			  	id      = elemNodes[i]
			  	mm[id]  += rho*N[i]*detJ*weights[gp]		
			  end
			end
	  	end
    end

    # time-independent Dirichlet boundary conditions on grid/solids
    fix_Dirichlet_grid(grid,data)
    fix_Dirichlet_solid(solids,data)


  #################################################
  # Time loop 
  #################################################

  while t < Tf

    @printf("Solving step: %d %f \n", counter, t)

    # ===========================================
    # reset grid data
    # ===========================================

    @inbounds for i = 1:nodeCount
	  nodalMass_S[i]      = 0.	  
	  nodalMomentum_S[i]  =  @SVector [0., 0.]
	  contact_solids[i]   = Set{Int64}()
    end

    boundary_nodes   = Set{Int64}()

    # ===========================================
    # particle to grid (deformable solids)
    # ===========================================

	@inbounds for s = 1:solidCount
		solid  = solids[s]
	
	    # deformable solids only
		if typeof(mats[s]) <: RigidMaterial continue end

		xx     = solid.pos
		mm     = solid.mass
		vv     = solid.velocity		
		fint   = solid.fint
		fbody  = solid.fbody
		du     = solid.dU

		nMass     = nodalMass[s]
		nMomenta0 = nodalMomentum0[s]		
		nMomenta  = nodalMomentum[s]		
		nForce    = nodalForce[s]
		normals   = nodalNormals[s]
		bnd_particles = solid.mesh.node_sets["boundary"]
		

        # reset body mass/momenta0/force 
		nMass[:]  .= 0.
		fill!(nMomenta0,@SVector[0.,0.])
		fill!(nForce,   @SVector[0.,0.])
		fill!(normals,  @SVector[0.,0.])
		

	  	@inbounds for ip = 1:solid.nodeCount
	        #getShapeAndGradient(nearPoints,funcs,ders,xx[ip], grid) 
			support   = getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid,basis)
	        
	        fMass     = mm[ip]
	        vp        = vv[ip]
	        fp        = fint[ip]
	        fb        = fbody[ip]
			
			#println(funcs)
			@inbounds for i = 1:support
				id    = nearPoints[i]; # index of node 'i'
				Ni    = funcs[i]				
				dNi   = @view ders[:,i]			
				Nim   = Ni * fMass
				# mass, momentum, internal force and external force
				nMass[id]      += Nim
				nMomenta0[id]  += Nim * vp 
				nForce[id]     -= Ni  * fp
				nForce[id]     += Ni  * fb
				normals[id]    += fMass * dNi
				if (ip in bnd_particles ) push!(boundary_nodes,id) end
				push!(contact_solids[id],s)  # add solid 's' to the set of contact solids at node 'id'
			end
			du[ip] = @SVector [0., 0.]
	  	end	  	
	
		# ===========================================
		# update grid
		# ===========================================

		@inbounds for i=1:nodeCount

		    nMomenta[i]          = nMomenta0[i] + nForce[i] * dtime
		    nodalMass_S[i]      += nMass[i]
		    nodalMomentum_S[i]  += nMomenta[i]

		    nMomenta0[i]  /= nMass[i]
		    nMomenta[i]   /= nMass[i]
			  
    
	        fixed_dirs       = @view grid.fixedNodes[:,i]
	        if fixed_dirs[1] == 1
	        	nMomenta0[i] = setindex(nMomenta0[i],0.,1)
	   		    nMomenta[i]  = setindex(nMomenta[i], 0.,1)
	        end
	        if fixed_dirs[2] == 1
				nMomenta0[i] = setindex(nMomenta0[i],0.,2)
				nMomenta[i]  = setindex(nMomenta[i], 0.,2)
	        end	         
		end
	end # end of loop over solids

	if haskey(data, "rigid_body_velo") 

	    rigid_solids = data["rigid_body_velo"]

	 	for (s,f) in rigid_solids
	 		solid       = solids[s]
			xx          = solid.centroids		
			surf_normals= solid.normals	
			normals     = nodalNormals[s]
	
			fill!(normals,  @SVector[0.,0.])
			vex,vey     = f(t)			
			
			@inbounds for ip = 1:solid.parCount
				#getAdjacentGridPoints(nearPointsLin,xx[ip],grid,linBasis)
				#support   = getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid,basis)
				support    = getShapeFunctions(nearPoints,funcs,xx[ip],grid, basis)
				#println(nearPoints)
				normal_ip = surf_normals[ip]
				@inbounds for i = 1:4
					id                  = nearPoints[i]; # index of node 'i'					
					#dNi                 = @view ders[:,i]		
					normals[id]        += funcs[i] * normal_ip
					mi                  = nodalMass_S[id]
					#if (ip in bnd_particles ) push!(boundary_nodes,id) end
				    push!(contact_solids[id],s)  # add solid 's' to the set of contact solids at node 'id'
	
				    nodalMomentum_S[id]   = setindex(nodalMomentum_S[id], mi*vex,1)
				    nodalMomentum_S[id]   = setindex(nodalMomentum_S[id], mi*vey,2)				    			    
				end
			end
			#solid.reaction_forces .= @SVector[Fx, Fy, Fz]
		end
	end

    ###########################################################################
    # contact treatments
    ###########################################################################
	
	# loop over all boundary nodes, some of them, which are contact nodes, need 
    # contact treatments
	@inbounds for i in boundary_nodes
	    solids_at_node_i = collect(contact_solids[i]) # to convert from set to array to access using index
	    
	    # skip nodes receive contribution from 1 single solid,
	    # as they are not contact nodes
	    if length(solids_at_node_i) == 1 continue end 

        # special treatment of normals here, average normal of solid 1 and solid 2
        # just for 2 deformable solids
        s1        = solids_at_node_i[1]
        s2        = solids_at_node_i[2]

#@printf("Contacting solids at node i: %d %d %d\n", s1, s2, i)

	    nI_body_1 = nodalNormals[s1][i]
	    nI_body_2 = nodalNormals[s2][i]

	    # if one solid is rigid, then use normal of it
	    # otherwise, use average

	    if     typeof(mats[s1]) <: RigidMaterial 
	    	nI_common .= nI_body_1 
	    	nI_common /= norm(nI_common)

		    #normals    = [nI_common -nI_common]
		    normals    = SMatrix{2,2}( nI_common[1], nI_common[2],
		    	                      -nI_common[1],-nI_common[2])
		    
#println(normals)
	    elseif typeof(mats[s2]) <: RigidMaterial 
	    	nI_common .= nI_body_2
	    	nI_common /= norm(nI_common)	    

		    #normals    = [nI_common -nI_common]
		    normals    = SMatrix{2,2}( -nI_common[1], -nI_common[2],
		    	                        nI_common[1], nI_common[2])	

	    else
	    	nI_common .= nI_body_1 .- nI_body_2
	    	nI_common /= norm(nI_common)

		    #normals    = [nI_common -nI_common]
		    normals    = SMatrix{2,2}( nI_common[1], nI_common[2],
		    	                      -nI_common[1],-nI_common[2])
	    end


	    for is = 1 : length(solids_at_node_i)
	    	s               = solids_at_node_i[is]
            # do not correct velocity of rigid body
	    	if typeof(mats[s]) <: RigidMaterial continue end

		    solid           = solids[s]
			nodalMomentum_1 = nodalMomentum[s]
			nodalMass_1     = nodalMass[s]
			

			velo1    = nodalMomentum_1[i];                       # body 1 velo				
		    velocm   = nodalMomentum_S[i]/nodalMass_S[i];        # system velo
		    
		    nI      .= normals[:,is]
		    deltaVe1 .= velo1 .- velocm;			    
		    D1       = dot(deltaVe1,  nI);			    

		    if ( D1 > 0 )
			    C1       =  deltaVe1[1]*nI[2] - deltaVe1[2]*nI[1];			    
			    absC1    = abs(C1);			    
			    muPrime1 = min(fric,absC1/D1);
		        
		        nodalMomentum_1[i] = velo1 - D1*(  nI + (muPrime1/absC1)*@SVector[ nI[2]*C1, - nI[1]*C1] );		        	        
		    
		        #println("approaching\n")
		    else			    	
		        #println("separating\n")
		    end		
	    end
	end

    # ====================================================================
    # grid to particle (deformable solids): update particle velocity/pos
    # ====================================================================

   @inbounds for s = 1:solidCount
		# only deformable solids here
	  	solid = solids[s]

	  	if typeof(mats[s]) <: RigidMaterial continue end

	  	xx    = solid.pos
	  	mm    = solid.mass
	  	vv    = solid.velocity
	  	du    = solid.dU
	    fix   = solid.fixedNodes
	  	nMomentum0 = nodalMomentum0[s]
	  	nMomentum  = nodalMomentum[s]
	  	nMass      = nodalMass[s]
	  	
	  	
	  	@inbounds for ip = 1:solid.nodeCount
			support   = getShapeFunctions(nearPoints,funcs,ip, grid, solid, basis)	        
			vvp       = vv[ip]
			xxp       = xx[ip]
			dup       = du[ip]
		
			for i = 1:support
				in = nearPoints[i]; # index of node 'i'
				Ni = funcs[i]
				mI = nMass[in]
			    if (mI > alg.tolerance)
					#invM       = 1.0 / mI
					vI         = nMomentum[in] #* invM
					#vvt        += Ni * vI  => too much dissipation
					vvp       += Ni * (nMomentum[in] - nMomentum0[in])# * invM
					xxp       += Ni * vI * dtime				
					dup       += Ni * vI * dtime				
		   		end
		   	end
			vv[ip]      = vvp
			xx[ip]      = xxp	
			du[ip]      = dup
		    # Dirichlet BCs on the mesh
			fixed_dirs       = @view fix[:,ip]
	        if fixed_dirs[1] == 1
				vv[ip] = setindex(vv[ip],0.,1)
				du[ip] = setindex(du[ip],0.,1)
				xx[ip] = setindex(xx[ip],solid.pos0[ip][1],1)
	        end
	        if fixed_dirs[2] == 1
				vv[ip] = setindex(vv[ip],0.,2)
				du[ip] = setindex(du[ip],0.,2)
				xx[ip] = setindex(xx[ip],solid.pos0[ip][2],2)
	        end
	    end
	end


    # ====================================================================
    #  update particle internal forces (FEM, so no contact here)
    # ====================================================================
t       += dtime
    @inbounds for s = 1:solidCount
		# only deformable solids here
	  	solid = solids[s]
	  	if typeof(mats[s]) <: RigidMaterial continue end

	  	xx     = solid.pos
	  	XX     = solid.pos0
	  	mm     = solid.mass
	  	du     = solid.dU
	  	velo   = solid.velocity
	  	F      = solid.deformationGradient
	  	mat    = mats[s]
	  	stress = solid.stress
	  	strain = solid.strain
	  	elems  = solid.elems
	  	fint   = solid.fint
	  	fbody  = solid.fbody
	  	vol    = solid.volume
	  	meshBasis = solid.basis

	  	for ip=1:solid.nodeCount 
	  		fint[ip]  = @SVector [0., 0.]
	  		fbody[ip] = @SVector [0., 0.]
	  	end
	  	g .= 0.
        # loop over solid elements, not solid nodes
	  	@inbounds for ip = 1:solid.parCount
			elemNodes =  @view elems[ip,:]  
			coords    =  @view XX[elemNodes]
	        vel_grad  =  SMatrix{2,2}(0., 0., 0., 0.)
	        vel_gradT =  SMatrix{2,2}(0., 0., 0., 0.)
			
			detJ      = lagrange_basis_derivatives!(N, dNdx, meshBasis, gpCoords, coords)
			w         = detJ * wgt
			vol[ip]   = w
			#println(dNdx)
			#println(sum(dNdx, dims=2))
			xcenter = ycenter = 0.
			for i = 1:length(elemNodes)
				in         = elemNodes[i]; # index of node 'i'
			    dNi        = @view dNdx[:,i]			
				vI         = du[in]
				vIt        = xx[in] - XX[in]
				#vI         = velo[in]
				vel_grad  += SMatrix{2,2}(dNi[1]*vI[1], dNi[1]*vI[2],
   										  dNi[2]*vI[1], dNi[2]*vI[2])
				vel_gradT += SMatrix{2,2}(dNi[1]*vIt[1], dNi[1]*vIt[2],
   										  dNi[2]*vIt[1], dNi[2]*vIt[2])

                xcenter   += N[i] * coords[i][1]			
                ycenter   += N[i] * coords[i][2]			
		   	end
			   	
		   	
		   	#F[ip]        += dtime * vel_grad
		   	# D            = 0.5 * (vel_grad + vel_grad')
		   	Fdot   .= vel_grad/dtime
		   	F[ip]      += Fdot*dtime
		   	#F[ip]        = Identity + vel_gradT
		   	#strain[ip]   = 0.5 * (vel_grad + vel_grad' + vel_grad * vel_grad')

		   	L      = Fdot*inv(F[ip])
		   	D      = 0.5 * dtime * (L + L')
		   	strain[ip]   += D #0.5 * (vel_grad + vel_grad' + vel_grad * vel_grad')


		   	J            = det(F[ip])
		   	if ( J < 0. )
					@printf("Troubled particle: %d %f %f \n", s, XX[ip][1], XX[ip][2])
					#println(F[ip])
					#@error("J is negative\n")
			end
		   	#println(strain[ip])
	   	     #@timeit "3" update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
	   	    update_stress!(stress[ip],mat,strain[ip],strain[ip] ,F[ip],J,ip,dtime)

            if s == 1 body(g,@SVector[0,0],t)  end
            sigma = stress[ip]
            P     = J*sigma*inv(F[ip])'  # convert to 1st Piola Kirchoof stress
            # compute nodal internal force fint
		    for i = 1:length(elemNodes)
				in  = elemNodes[i]; # index of node 'i'
			    dNi = @view dNdx[:,i]			
	   	        fint[in]     += w * @SVector[P[1,1] * dNi[1] + P[1,2] * dNi[2],
											 P[2,1] * dNi[1] + P[2,2] * dNi[2]]
                #fbody[in]    += w*mat.density*N[i]*g# because N_1=N_2==.25													        
                fbody[in]    += w*mat.density*N[i]*g															        
            end
	   end
	end

	##################################################
	# update position of rigid solids
    ##################################################
	
	if haskey(data, "rigid_body_velo") 

	    rigid_solids = data["rigid_body_velo"]

	 	for (s,f) in rigid_solids
	 		solid       = solids[s]
			xx          = solid.pos
			xc          = solid.centroids
			(vex,vey) = f(t)

			if (vex,vey) == (0.,0.) continue end
			
			@inbounds for ip = 1:solid.nodeCount
		      xx[ip]   += dtime * @SVector [vex,vey]
		    end

		    @inbounds for ip = 1:solid.parCount
		      xc[ip]   += dtime * @SVector [vex,vey]
		    end
		end
	end


	##################################################
	# output
    ##################################################

	if (counter%output.interval == 0)
		plotParticles_2D(output,solids,mats,counter)
		compute_femp(fixes,t)
	end

#    t       += dtime
    counter += 1
  end # end of time loop
  #@printf("Solving done \n")
  closeFile(fixes)
end # end solve()




