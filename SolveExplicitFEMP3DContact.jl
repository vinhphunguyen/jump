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
function solve_explicit_dynamics_femp_3D_Contact(grid,solids,mats,basis,body,fric,alg::TLFEM,output,fixes,Tf,dtime)
    t       = 0.
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

    boundary_nodes   = Vector{Set{Int64}}(undef,solidCount)



    # pre_allocating arrays for temporary variable
   
	nearPoints,funcs, ders = initialise(grid,basis)

    # for rigid bodies
	linBasis       = LinearBasis()
	nearPointsLin  = [0,0,0,0]

    # temporary variables, just 1 memory allocation
    Fdot      = zeros(3,3)
    vel_grad  = SMatrix{3,3}(0,0,0,0,0,0,0,0,0)
    D         = SMatrix{3,3}(0,0,0,0,0,0,0,0,0) #zeros(Float64,2,2)
    g         = [0.,0.,0.]
    deltaVe1  = zeros(3)
    dVxNI     = zeros(3)
    omega     = zeros(3)
    nI        = zeros(3)
   

    @inbounds for s = 1:solidCount
		solid  = solids[s]
        initializeBasis(solid,mats[s].density)
    end

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
    end

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
		bnd_nodes     = Set()#solid.boundary_nodes
		

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
				if (ip in bnd_particles ) push!(bnd_nodes,id) end
			end
			du[ip] = @SVector [0., 0., 0.]
	  	end

	  	boundary_nodes[s] = bnd_nodes
	
		# ===========================================
		# update grid
		# ===========================================

		@inbounds for i=1:nodeCount

		    nMomenta[i]          = nMomenta0[i] + nForce[i] * dtime
		    nodalMass_S[i]      += nMass[i]
		    nodalMomentum_S[i]  += nMomenta[i]

		    nMomenta0[i]  /= nMass[i]
		    nMomenta[i]   /= nMass[i]
			  

	        # apply Dirichet boundary conditions on which solid?
	        
			if grid.fixedXNodes[i] == 1
				nMomenta0[i] = setindex(nMomenta0[i],0.,1)
	   		    nMomenta[i]  = setindex(nMomenta[i], 0.,1)
			end
	        if grid.fixedYNodes[i] == 1
				nMomenta0[i] = setindex(nMomenta0[i],0.,2)
				nMomenta[i]  = setindex(nMomenta[i], 0.,2)
	        end	 
	        if grid.fixedZNodes[i] == 1
				nMomenta0[i] = setindex(nMomenta0[i],0.,3)
				nMomenta[i]  = setindex(nMomenta[i], 0.,3)
	        end	 
	         
		end
	end # end of loop over solids

	# ===========================================
	# particle to grid (rigid solids)
	# ===========================================

# 	@inbounds for s = 1:solidCount
# 		solid  = solids[s]
# 		# rigid solids only
# 		if !(typeof(mats[s]) <: RigidMaterial) continue end
# 		xx      = solid.pos
# 		vex     = solid.mat.vx
# 		vey     = solid.mat.vy
# 		@inbounds for ip = 1:solid.parCount
# 			getAdjacentGridPoints(nearPointsLin,xx[ip],grid,linBasis)
# #			println(nearPoints)
# 			@inbounds for i = 1:4
# 				id                  = nearPointsLin[i]; # index of node 'i'
# 				mi                  = nodalMass[id]
# 				if solid.mat.fixed
# 				  nodalMomentum[id]   = setindex(nodalMomentum[id],0.,1)
# 				  nodalMomentum[id]   = setindex(nodalMomentum[id],0.,2)
#   				  nodalMomentum0[id]  = setindex(nodalMomentum0[id],0.,1)
#   				  nodalMomentum0[id]  = setindex(nodalMomentum0[id],0.,2)
# 			    else
# 					if vex != 0.
# 					  nodalMomentum[id]   = setindex(nodalMomentum[id], mi*vex,1)
# 					  nodalMomentum0[id]  = setindex(nodalMomentum0[id],mi*vex,1)
# 				    end
# 					if vey != 0.
# 					  nodalMomentum[id]   = setindex(nodalMomentum[id], mi*vey,2)
# 					  nodalMomentum0[id]  = setindex(nodalMomentum0[id],mi*vey,2)
# 					end
# 				end
# 				#println(nodalMomentum)
# 			end
# 		end
# 	end

    ###########################################################################
    # contact treatments
    ###########################################################################
	@inbounds for s = 1:solidCount
		solid           = solids[s]
		nodalMomentum_1 = nodalMomentum[s]
		nodalMass_1     = nodalMass[s]
		bnd_nodes       = boundary_nodes[s]
		normals         = nodalNormals[s]
		@inbounds for i in bnd_nodes
			if ( nodalMass_S[i] != nodalMass_1[i] ) # this is  a contact node
				println(i)
				velo1    = nodalMomentum_1[i];        # body 1 velo				
			    velocm   = nodalMomentum_S[i]/nodalMass_S[i];        # system velo
			    nI      .= normals[i] / norm(normals[i])
			    deltaVe1 .= velo1 .- velocm;			    
			    D1       = dot(deltaVe1,  nI);			    
    
			    if ( D1 > 0 )
				    dVxNI   .= cross(deltaVe1,  nI);			    
				    C1       = norm(dVxNI)		    
				    omega   .= dVxNI / C1
				    muPrime1 = min(fric,C1/D1);
			        
			        nodalMomentum_1[i] = velo1 - D1*(  nI + muPrime1*cross(nI,omega) );		        
			    
			        println("approaching\n")
			    else			    	
			        println("separating\n")
			    end
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
	    fixX  = solid.fixedNodesX
	  	fixY  = solid.fixedNodesY
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
    #  update particle internal forces (FEM, so no contact here)
    # ====================================================================
t       += dtime
    
     @inbounds for s = 1:solidCount
		# only deformable solids here
	  	solid = solids[s]

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
				vIt        = xx[in] - XX[in]
				vel_grad  += SMatrix{3,3}(dNi[1]*vI[1], dNi[1]*vI[2], dNi[1]*vI[3],
   										  dNi[2]*vI[1], dNi[2]*vI[2], dNi[2]*vI[3],
   										  dNi[3]*vI[1], dNi[3]*vI[2], dNi[3]*vI[3] )

				vel_gradT += SMatrix{3,3}(dNi[1]*vIt[1], dNi[1]*vIt[2], dNi[1]*vIt[3],
   										  dNi[2]*vIt[1], dNi[2]*vIt[2], dNi[2]*vIt[3],
   										  dNi[3]*vIt[1], dNi[3]*vIt[2], dNi[3]*vIt[3] )
		   	end
		   	
			#dstrain      = 0.5 * (vel_grad + vel_grad' + vel_grad * vel_grad') - strain[ip] 
			
			Fdot         .= vel_grad/dtime		   
		   	F[ip]         = Identity + vel_gradT
		   	J             = det(F[ip])
		   	L             = Fdot*inv(F[ip])
		   	D             = 0.5 * dtime * (L + L')
		   	strain[ip]   += D #0.5 * (vel_grad + vel_grad' + vel_grad * vel_grad')

		   	#println(strain[ip])
	   	     #@timeit "3" update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
	   	    update_stress!(stress[ip],mat,strain[ip],D,F[ip],J,ip,dtime)

            sigma = stress[ip]
            P     = J*sigma*inv(F[ip])'  # convert to Piola Kirchoof stress

            body(g,[0 0 0],t)  
            # compute nodal internal force fint
		    for i = 1:length(elemNodes)
				in  = elemNodes[i]; # index of node 'i'
			    dNi = @view dNdx[:,i]			
	   	        fint[in]  +=  detJ * @SVector[P[1,1] * dNi[1] + P[1,2] * dNi[2] + P[1,3] * dNi[3],
										      P[2,1] * dNi[1] + P[2,2] * dNi[2] + P[2,3] * dNi[3],
										      P[3,1] * dNi[1] + P[3,2] * dNi[2] + P[3,3] * dNi[3] ]
                fbody[in] += detJ*mat.density*N[i]*g								        
            end
	   end
	end

	##################################################
	# update position of rigid solids
    ##################################################
	@inbounds for s = 1:solidCount
		# only rigid solids here
	  	solid = solids[s]
		if !(typeof(mats[s]) <: RigidMaterial) continue end

	  	xx    = solid.pos
	  	vx    = solid.mat.vx
	  	vy    = solid.mat.vy
		#println(ve)
        @inbounds for ip = 1:solid.nodeCount
	      xx[ip]   += dtime * @SVector [vx,vy]
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




