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
function solve_explicit_dynamics_femp_3D_Contact(grid,solids,mats,basis,body,alg::TLFEM,output,fixes,data)
    Tf    = data["total_time"]:: Float64
	dtime = data["dt"]        :: Float64    
	t     = data["time"]      :: Float64     
	fric  = data["friction"]  :: Float64

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

    # boundary nodes (nodes surrounding the surfaces of all solids)
    # a sub-set of them are contact nodes that requires contact treatment
    # old implementation: store for each solid
    #boundary_nodes   = Vector{Set{Int64}}(undef,solidCount)
    # new implementation:   
    #boundary_nodes   = Set{Int64}()

    # for all grid nodes, store ids of solids that are in contact
    contact_solids   = Vector{Set{Int64}}(undef,nodeCount)

    # pre_allocating arrays for temporary variable
   
	nearPoints,funcs, ders = initialise(grid,basis)

    # for rigid bodies
	linBasis       = LinearBasis()
	nearPointsLin  = [0,0,0,0,0,0,0,0]

    # temporary variables, just 1 memory allocation
    Fdot      = SMatrix{3,3}(0,0,0,0,0,0,0,0,0) #zeros(3,3)
    sigma     = SMatrix{3,3}(0,0,0,0,0,0,0,0,0) #zeros(3,3)
    vel_grad  = SMatrix{3,3}(0,0,0,0,0,0,0,0,0)
    D         = SMatrix{3,3}(0,0,0,0,0,0,0,0,0) #zeros(Float64,2,2)
    g         = [0.,0.,0.]
    deltaVe1  = zeros(3)
    dVxNI     = zeros(3)
    omega     = zeros(3)
    nI        = zeros(3)
    nI_common = zeros(3)
   
    # compute nodal mass and FE shape functions and derivatives at centroids for all
    # elements in the mesh of each solid
    @inbounds for s = 1:solidCount
		solid  = solids[s]
		if typeof(mats[s]) <: RigidMaterial 
			compute_normals!(solid)
			continue 
		end
        initializeBasis(solid,mats[s].density)
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
	  nodalMomentum_S[i]  =  @SVector [0., 0., 0.]
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
		fill!(nMomenta0,@SVector[0.,0.,0.])
		fill!(nForce,   @SVector[0.,0.,0.])
		fill!(normals,  @SVector[0.,0.,0.])
		#bnd_nodes     = Set{Int64}() # old implementation: solid.boundary_nodes
		
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
                nMass[id]      += Nim             # check: no memory alloc
                nMomenta0[id]  .+= Nim.* vp 
				nForce[id]     -= Ni  * fp
                nForce[id]    .+= Ni .* fb       # check: no memory alloc
                normals[id]    .+= fMass .* dNi
				#if (ip in bnd_particles ) push!(bnd_nodes,id) end
				if (ip in bnd_particles ) push!(boundary_nodes,id) end
				push!(contact_solids[id],s)  # add solid 's' to the set of contact solids at node 'id'
			end
			du[ip] = @SVector [0., 0., 0.]
	  	end

	  	#boundary_nodes[s] = bnd_nodes old implementation
	
		# ===========================================
		# update grid
		# ===========================================

		@inbounds for i=1:nodeCount

            nMomenta[i]         .= nMomenta0[i] .+ dtime .* nForce[i] 
            nodalMass_S[i]      += nMass[i]
		    nodalMomentum_S[i]  += nMomenta[i]

            nMomenta0[i]  ./= nMass[i]
		    nMomenta[i]   ./= nMass[i]
			  

	        # apply Dirichet boundary conditions on which solid?

	        
	        fixed_dirs       = @view grid.fixedNodes[:,i]
	        if fixed_dirs[1] == 1
                 #nMomenta0[i] = setindex(nMomenta0[i],0.,1) => memory alloc!!!
                nMomenta0[i][1] = 0.
                nMomenta[i][1]  = 0.
	   		    
	        end
	        if fixed_dirs[2] == 1
			    nMomenta0[i][2] = 0.
                nMomenta[i][2]  = 0.
	        end
	        if fixed_dirs[3] == 1
			    nMomenta0[i][3] = 0.
                nMomenta[i][3]  = 0.
	        end
	   
	         
		end
	end # end of loop over solids

	# ===========================================
	# particle to grid (rigid solids)
	# ===========================================

	if haskey(data, "rigid_body_velo") 

	    rigid_solids = data["rigid_body_velo"]::Array{Tuple{Int64,Function},1}

	 	for (s,f) in rigid_solids
	 		solid       = solids[s]
			xx          = solid.centroids		
			surf_normals= solid.normals	
			normals     = nodalNormals[s]
	
			fill!(normals,  @SVector[0.,0.,0.])
			(vex,vey,vez) = f(t)::Tuple{Float64, Float64, Float64} 

			#Fx = Fy = Fz = 0.
			# loop over centroids of surface elements
			@inbounds for ip = 1:solid.surfCount
			    support    = getShapeFunctions(nearPoints,funcs,xx[ip],grid, basis)
				#println(nearPoints)
				normal_ip = surf_normals[ip]
	#			println(nearPoints)
				@inbounds for i = 1:8
					id                  = nearPoints[i]; # index of node 'i'					
                    normals[id]       .+= funcs[i] .* normal_ip
					mi                  = nodalMass_S[id]
					#if (ip in bnd_particles ) push!(boundary_nodes,id) end
				    push!(contact_solids[id],s)  # add solid 's' to the set of contact solids at node 'id'

                    # Fx += (mi*vex - nodalMomentum_S[id][1])/dtime
                    # Fy += (mi*vey - nodalMomentum_S[id][2])/dtime
                    # Fz += (mi*vez - nodalMomentum_S[id][3])/dtime
	
                    nodalMomentum_S[id]   .=  mi*@SVector[vex,vey,vez]				    		    
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

        #@printf("Contacting solids: %d %d\n", s1, s2)

	    nI_body_1 = nodalNormals[s1][i]
	    nI_body_2 = nodalNormals[s2][i]
	    # if one solid is rigid, then use normal of it
	    # otherwise, use average

	    if     typeof(mats[s1]) <: RigidMaterial 
	    	nI_common .= nI_body_1 
            nI_common /= norm(nI_common)

		    #normals    = [nI_common -nI_common]
 		    normals    = SMatrix{3,2}( nI_common[1], nI_common[2], nI_common[3],
		    	                      -nI_common[1],-nI_common[2],-nI_common[3])

	    elseif typeof(mats[s2]) <: RigidMaterial 
	    	nI_common .= nI_body_2
	    	nI_common /= norm(nI_common)	    

		    #normals    = [nI_common -nI_common]
		    normals    = SMatrix{3,2}( -nI_common[1], -nI_common[2], -nI_common[3],
		    	                        nI_common[1], nI_common[2],   nI_common[3])	
	    else
	    	nI_common .= nI_body_1 .- nI_body_2
	    	nI_common /= norm(nI_common)

		    #normals    = [nI_common -nI_common]
		    normals    = SMatrix{3,2}( nI_common[1], nI_common[2], nI_common[3],
		    	                      -nI_common[1],-nI_common[2],-nI_common[3])
	    end


	    for is = 1 : length(solids_at_node_i)
	    	s               = solids_at_node_i[is]
            # do not correct velocity of rigid body
	    	if typeof(mats[s]) <: RigidMaterial continue end

		    solid           = solids[s]
			nodalMomentum_1 = nodalMomentum[s]
			nodalMass_1     = nodalMass[s]
			

			velo1    = nodalMomentum_1[i];                       # body 1 velo				
            velocm   = nodalMomentum_S[i]./nodalMass_S[i];        # system velo
		    
		    nI      .= normals[:,is]
		    deltaVe1 .= velo1 .- velocm;			    
		    D1       = dot(deltaVe1,  nI);			    

		    if ( D1 > 0 )
			    dVxNI   .= cross(deltaVe1,  nI);			    
			    C1       = norm(dVxNI)		    
			    omega   .= dVxNI / C1
			    muPrime1 = min(fric,C1/D1);
		        
		        nodalMomentum_1[i] = velo1 - D1*(  nI + muPrime1*cross(nI,omega) );		        
		    
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

	  	xx         = solid.pos
	  	mm         = solid.mass
	  	vv         = solid.velocity
	  	du         = solid.dU
	    fix        = solid.fixedNodes
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
	        if fixed_dirs[3] == 1
				vv[ip] = setindex(vv[ip],0.,3)
				du[ip] = setindex(du[ip],0.,3)
				xx[ip] = setindex(xx[ip],solid.pos0[ip][3],3)
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

	  	XX     = solid.pos0
	  	xx     = solid.pos
	  	mm     = solid.mass
	  	du     = solid.dU
	  	F      = solid.deformationGradient
	  	mat    = mats[s]
	  	stress = solid.stress
	  	strain = solid.strain
	  	elems  = solid.elems
	  	fint   = solid.fint
	  	fbody  = solid.fbody
	  	vol    = solid.volume
	  	meshBasis = solid.basis
	  	feFuncs   = solid.N
	  	feDers    = solid.dNdx
	  	feJaco    = solid.detJ
	  	
	  	for ip=1:solid.nodeCount 
	  		fint[ip]  = @SVector [0., 0., 0.]
	  		fbody[ip] = @SVector [0., 0., 0.]
	  	end
        # loop over solid elements, not solid nodes
	  	@inbounds for ip = 1:solid.parCount
			elemNodes =  @view elems[ip,:]  
			#coords    =  @view XX[elemNodes]
	        vel_grad   =  SMatrix{3,3}(0,0,0,0,0,0,0,0,0)
	        vel_gradT  =  SMatrix{3,3}(0,0,0,0,0,0,0,0,0)
			
			#detJ      = lagrange_basis_derivatives!(N, dNdx, meshBasis, gpCoords, coords)
			detJ      = feJaco[ip]
			N         = feFuncs[ip]
			dNdx      = feDers[ip]
			
			vol[ip]   = detJ
			#println(dNdx)
			#println(sum(dNdx, dims=2))
			for i = 1:length(elemNodes)
				in         = elemNodes[i]; # index of node 'i'
			    dNi        = @view dNdx[:,i]			
				vI         = du[in]
                vIt        = xx[in] - XX[in]   # no memory alloc
				vel_grad  += SMatrix{3,3}(dNi[1]*vI[1], dNi[1]*vI[2], dNi[1]*vI[3],
   										  dNi[2]*vI[1], dNi[2]*vI[2], dNi[2]*vI[3],
   										  dNi[3]*vI[1], dNi[3]*vI[2], dNi[3]*vI[3] )

				vel_gradT += SMatrix{3,3}(dNi[1]*vIt[1], dNi[1]*vIt[2], dNi[1]*vIt[3],
   										  dNi[2]*vIt[1], dNi[2]*vIt[2], dNi[2]*vIt[3],
   										  dNi[3]*vIt[1], dNi[3]*vIt[2], dNi[3]*vIt[3] )
		   	end
		   	
			#dstrain      = 0.5 * (vel_grad + vel_grad' + vel_grad * vel_grad') - strain[ip] 
			
			Fdot          = vel_grad/dtime		   
       	    F[ip]         = Identity + vel_gradT  # no memory alloc
		   	J             = det(F[ip])
            Finv          = inv(F[ip])    # no memory alloc
            L             = Fdot*Finv     # no memory alloc when Fdot is a SMatrix
#@timeit "2"            D             .= (0.5 * dtime) .* (L .+ L')
            D             = (0.5 * dtime) * (L + L')  # no memory alloc
		   	strain[ip]   += D #0.5 * (vel_grad + vel_grad' + vel_grad * vel_grad')

		   	#println(strain[ip])
	   	     #@timeit "3" update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
	   	    update_stress!(stress[ip],mat,strain[ip],D,F[ip],J,ip,dtime)

            #sigma = SMatrix{3,3}(0,0,0,0,0,0,0,0,0) #stress[ip]
            sigma = stress[ip]
            P     = J*sigma*Finv'  # convert to Piola Kirchoof stress, no memory alloc

            #body(g,[0 0 0],t)  
            # compute nodal internal force fint
		    for i = 1:length(elemNodes)
				in  = elemNodes[i]; # index of node 'i'
			    dNi = @view dNdx[:,i]			
                fint[in]  +=  detJ * @SVector[P[1,1] * dNi[1] + P[1,2] * dNi[2] + P[1,3] * dNi[3],
										      P[2,1] * dNi[1] + P[2,2] * dNi[2] + P[2,3] * dNi[3],
										      P[3,1] * dNi[1] + P[3,2] * dNi[2] + P[3,3] * dNi[3] ]
# Performance issue here: mat.density is ANY because we have heterogeneous contains of mats: RIgid/other mats.										      
#@timeit "4"                fbody[in] += detJ*mat.density*N[i]*g								        
            end
	   end
	end

	fix_Dirichlet_solid(solids,data,dtime)

	##################################################
	# update position of rigid solids
    ##################################################
	
		
	if haskey(data, "rigid_body_velo") 
	    
	    rigid_solids = data["rigid_body_velo"]::Array{Tuple{Int64,Function},1}

	 	for (s,f) in rigid_solids
	 		solid       = solids[s]
			xx          = solid.pos
			xc          = solid.centroids
			(vex,vey,vez) = f(t)::Tuple{Float64, Float64, Float64} 

			if (vex,vey,vez) == (0.,0.,0.) continue end
			# update all nodes for visualization only
			@inbounds for ip = 1:solid.nodeCount
		      xx[ip]   += dtime * @SVector [vex,vey,vez]   # no memory alloc, 
		    end
            # update the centroids of surface elems
		    @inbounds for ip = 1:solid.surfCount
		      xc[ip]   += dtime * @SVector [vex,vey,vez]
		    end
		end
	end

	if (counter%output.interval == 0)
		plotParticles_3D(output,solids,mats,counter)
		compute_femp(fixes,t)
	end

#    t       += dtime
    counter += 1
  end # end of time loop
  #@printf("Solving done \n")
  closeFile(fixes)
end # end solve()




