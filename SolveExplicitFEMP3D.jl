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
function solve_explicit_dynamics_femp_3D(grid,solids,basis,alg::USL,output,fixes,Tf,dtime)
    t       = 0.
    counter = 0

    Identity       = UniformScaling(1.)
	solidCount     = length(solids)
	nodalMass      = grid.mass
	nodalMomentum0 = grid.momentum0
	nodalMomentum  = grid.momentum
	nodalForce     = grid.force
	g              = 1.#-1e6#grid.gravity

    D              = SMatrix{3,3}(0,0,0,0,0,0,0,0,0) #zeros(Float64,2,2)

    # pre_allocating arrays for temporary variable

	nearPoints,funcs, ders = initialise(grid,basis)

	nodePerElem = size(solids[1].elems,2)

	if nodePerElem == 4 
		meshBasis = Tet4()  
		wgt       = 1.
		weights   = [1.]
		gpCoords  = [0.25,0.25,0.25]
	end
	if nodePerElem == 8 
		meshBasis = Hexa8() 
		wgt       = 8.0
        weights   = [8.]
        gpCoords  = [0.,0.,0.]
	end

    noGP      = 1
    dNdx      = zeros(3,nodePerElem)
    N         = zeros(nodePerElem)#@SVector [0,0,0,0]
    vel_grad  = SMatrix{3,3}(0,0,0,0,0,0,0,0,0)

    # compute nodal mass (only once)

    # gpCoords  = zeros(3,8)
    # weights   = ones(8)
    # gpCoords[1,1] = -0.5773502691896257;gpCoords[2,1] = -0.5773502691896257;gpCoords[3,1] = -0.5773502691896257
    # gpCoords[1,2] =  0.5773502691896257;gpCoords[2,2] = -0.5773502691896257;gpCoords[3,2] = -0.5773502691896257
    # gpCoords[1,3] =  0.5773502691896257;gpCoords[2,3] =  0.5773502691896257;gpCoords[3,3] = -0.5773502691896257
    # gpCoords[1,4] = -0.5773502691896257;gpCoords[2,4] =  0.5773502691896257;gpCoords[3,4] = -0.5773502691896257
    # gpCoords[1,5] = -0.5773502691896257;gpCoords[2,5] =  0.5773502691896257;gpCoords[3,5] =  0.5773502691896257
    # gpCoords[1,6] = -0.5773502691896257;gpCoords[2,6] =  0.5773502691896257;gpCoords[3,6] =  0.5773502691896257
    # gpCoords[1,7] = -0.5773502691896257;gpCoords[2,7] =  0.5773502691896257;gpCoords[3,7] =  0.5773502691896257
    # gpCoords[1,8] = -0.5773502691896257;gpCoords[2,8] =  0.5773502691896257;gpCoords[3,8] =  0.5773502691896257

    # noGP = 8

    @inbounds for s = 1:solidCount
		solid  = solids[s]
		xx     = solid.pos
		mm     = solid.mass  # to be updated here
		elems  = solid.elems
		mat    = solid.mat
		vol    = solid.volume

	  	@inbounds for ip = 1:solid.parCount
			elemNodes  =  @view elems[ip,:]  
			elemNodes0 =   elems[ip,:]  
			coords     =  @view xx[elemNodes]
	    
			@inbounds for gp = 1:noGP
			  xieta = @view gpCoords[:,gp]
			  detJ  = lagrange_basis!(N, meshBasis, xieta, coords)
			  if detJ < 0.
			  	@printf("Negative Jacobian in Gmsh mesh file!!! \n")
			  	elemNodes[2] = elemNodes0[4]			  	
			  	elemNodes[4] = elemNodes0[2]
			  end
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
	  nodalMomentum0[i] =  @SVector [0., 0., 0.]
	  nodalForce[i]     =  @SVector [0., 0., 0.]
    end

    # ===========================================
    # particle (nodes) to grid (deformable solids)
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
			du[ip] = @SVector [0., 0., 0.]
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
    #  update particle internal forces ()
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
	  	fbody  = solid.fbody
	  	vol    = solid.volume
	  	for ip=1:solid.nodeCount 
	  		fint[ip]  = @SVector [0., 0., 0.]
	  		fbody[ip] = @SVector [0., 0., 0.]
	  	end
        # loop over solid elements, not solid nodes
	  	@inbounds for ip = 1:solid.parCount
			elemNodes =  @view elems[ip,:]  
			coords    =  @view xx[elemNodes]
	        vel_grad  =  SMatrix{3,3}(0,0,0,0,0,0,0,0,0)
			
			detJ      = lagrange_basis_derivatives!(dNdx, meshBasis, gpCoords, coords)
			vol[ip]   = detJ * 8 # weight of GP = 8
			#println(dNdx)
			#println(sum(dNdx, dims=2))
			for i = 1:length(elemNodes)
				in         = elemNodes[i]; # index of node 'i'
			    dNi        = @view dNdx[:,i]			
				vI         = du[in]
				vel_grad  += SMatrix{3,3}(dNi[1]*vI[1], dNi[1]*vI[2], dNi[1]*vI[3],
   										  dNi[2]*vI[1], dNi[2]*vI[2], dNi[2]*vI[3],
   										  dNi[3]*vI[1], dNi[3]*vI[2], dNi[3]*vI[3] )
		   	end
	
		   	D           = 0.5 * (vel_grad + vel_grad')
		   	strain[ip]  += D
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
	   	        fint[in]+=detJ*wgt*@SVector[sigma[1,1] * dNi[1] + sigma[1,2] * dNi[2] + sigma[1,3] * dNi[3],
										    sigma[1,2] * dNi[1] + sigma[2,2] * dNi[2] + sigma[2,3] * dNi[3],
										    sigma[1,3] * dNi[1] + sigma[2,3] * dNi[2] + sigma[3,3] * dNi[3] ]
                fbody[in] += detJ*mat.density*@SVector[0.,g,0.]										  
            end
	   end
	end

	if (counter%output.interval == 0)
		plotParticles_3D(output,solids,counter)
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
function solve_explicit_dynamics_femp_3D(grid,solids,basis,alg::TLFEM,output,fixes,Tf,dtime)
    t       = 0.
    counter = 0

    Identity       = UniformScaling(1.)
	solidCount     = length(solids)
	nodalMass      = grid.mass
	nodalMomentum0 = grid.momentum0
	nodalMomentum  = grid.momentum
	nodalForce     = grid.force

	g              = -1e6

    D              = SMatrix{3,3}(0,0,0,0,0,0,0,0,0) #zeros(Float64,2,2)

    # pre_allocating arrays for temporary variable
   
	nearPoints,funcs, ders = initialise(grid,basis)

    nodePerElem = size(solids[1].elems,2)

	if nodePerElem == 4 
		meshBasis = Tet4()  
		wgt       = 1.
		weights   = [1.]
		gpCoords  = [0.25,0.25,0.25]
	end
	if nodePerElem == 8 
		meshBasis = Hexa8() 
		wgt       = 8.0
        weights   = [8.]
        gpCoords  = [0.,0.,0.]
	end

    noGP      = 1
    dNdx      = zeros(3,nodePerElem)
    N         = zeros(nodePerElem)#@SVector [0,0,0,0]

    vel_grad  = SMatrix{3,3}(0,0,0,0,0,0,0,0,0)
   

    # compute nodal mass (only once)

    # gpCoords  = zeros(3,8)
    # weights   = ones(8)
    # gpCoords[1,1] = -0.5773502691896257;gpCoords[2,1] = -0.5773502691896257;gpCoords[3,1] = -0.5773502691896257
    # gpCoords[1,2] =  0.5773502691896257;gpCoords[2,2] = -0.5773502691896257;gpCoords[3,2] = -0.5773502691896257
    # gpCoords[1,3] =  0.5773502691896257;gpCoords[2,3] =  0.5773502691896257;gpCoords[3,3] = -0.5773502691896257
    # gpCoords[1,4] = -0.5773502691896257;gpCoords[2,4] =  0.5773502691896257;gpCoords[3,4] = -0.5773502691896257
    # gpCoords[1,5] = -0.5773502691896257;gpCoords[2,5] =  0.5773502691896257;gpCoords[3,5] =  0.5773502691896257
    # gpCoords[1,6] = -0.5773502691896257;gpCoords[2,6] =  0.5773502691896257;gpCoords[3,6] =  0.5773502691896257
    # gpCoords[1,7] = -0.5773502691896257;gpCoords[2,7] =  0.5773502691896257;gpCoords[3,7] =  0.5773502691896257
    # gpCoords[1,8] = -0.5773502691896257;gpCoords[2,8] =  0.5773502691896257;gpCoords[3,8] =  0.5773502691896257


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
			  detJ  = lagrange_basis!(N, meshBasis, xieta, coords)
			  if detJ < 0.
			  	elemNodes[2] = elemNodes0[4]			  	
			  	elemNodes[4] = elemNodes0[2]
			  end
			  #println(N)	
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
	  nodalMomentum0[i] =  @SVector [0., 0., 0.]
	  nodalForce[i]     =  @SVector [0., 0., 0.]
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
			support   = getShapeFunctions(nearPoints,funcs,ip, grid, solid,basis)
	        
	        fMass     = mm[ip]
	        vp        = vv[ip]
	        fp        = fint[ip]
	        fb        = fbody[ip]
			#body      = problem.bodyforce(xx[ip],t)
			
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
        if grid.fixedZNodes[i] == 1
			nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,3)
			nodalMomentum[i]  = setindex(nodalMomentum[i], 0.,3)
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
	  	fixZ  = solid.fixedNodesZ
	  	
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

    @inbounds for s = 1:solidCount
		# only deformable solids here
	  	solid = solids[s]

	  	XX     = solid.pos0
	  	mm     = solid.mass
	  	du     = solid.dU
	  	F      = solid.deformationGradient
	  	mat    = solid.mat
	  	stress = solid.stress
	  	strain = solid.strain
	  	elems  = solid.elems
	  	fint   = solid.fint
	  	fbody  = solid.fbody
	  	vol    = solid.volume
	  	for ip=1:solid.nodeCount 
	  		fint[ip]  = @SVector [0., 0., 0.]
	  		fbody[ip] = @SVector [0., 0., 0.]
	  	end
        # loop over solid elements, not solid nodes
	  	@inbounds for ip = 1:solid.parCount
			elemNodes =  @view elems[ip,:]  
			coords    =  @view XX[elemNodes]
	        vel_grad  =  SMatrix{3,3}(0,0,0,0,0,0,0,0,0)
			
			detJ      = lagrange_basis_derivatives!(N, dNdx, meshBasis, gpCoords, coords)
			w         = detJ * wgt
			vol[ip]   = w
			#println(dNdx)
			#println(sum(dNdx, dims=2))
			for i = 1:length(elemNodes)
				in         = elemNodes[i]; # index of node 'i'
			    dNi        = @view dNdx[:,i]			
				vI         = du[in]
				vel_grad  += SMatrix{3,3}(dNi[1]*vI[1], dNi[1]*vI[2], dNi[1]*vI[3],
   										  dNi[2]*vI[1], dNi[2]*vI[2], dNi[2]*vI[3],
   										  dNi[3]*vI[1], dNi[3]*vI[2], dNi[3]*vI[3] )
		   	end
	
		   	D            = 0.5 * (vel_grad + vel_grad')
		   	strain[ip]   = D
		   	F[ip]        = Identity + vel_grad
		   	J            = det(F[ip])
		   	#println(strain[ip])
	   	     #@timeit "3" update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
	   	    update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)

            sigma = stress[ip]
            P     = J*sigma*inv(F[ip])'  # convert to Piola Kirchoof stress

            # compute nodal internal force fint
		    for i = 1:length(elemNodes)
				in  = elemNodes[i]; # index of node 'i'
			    dNi = @view dNdx[:,i]			
	   	        fint[in]  +=  w * @SVector[P[1,1] * dNi[1] + P[1,2] * dNi[2] + P[1,3] * dNi[3],
										   P[2,1] * dNi[1] + P[2,2] * dNi[2] + P[2,3] * dNi[3],
										   P[3,1] * dNi[1] + P[3,2] * dNi[2] + P[3,3] * dNi[3] ]
                fbody[in] += w*mat.density*N[i]*@SVector[0.,g,0.]								        
            end
	   end
	end

	if (counter%output.interval == 0)
		plotParticles_3D(output,solids,counter)
		plotGrid(output,grid,counter)
		compute_femp(fixes,t)
	end

    t       += dtime
    counter += 1
  end # end of time loop
  #@printf("Solving done \n")
  closeFile(fixes)
end # end solve()
