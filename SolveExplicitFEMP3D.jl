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

# This file contains functions for solving 3D explicit dynamics problems
# using the Generalized Particle in Cell (GPIC) method
# 29 April 2022: 
#  --  ULFEM, TLFEM with one point quad and TLFEMFull with full quadrature
#  -- AxiSymmetric




######################################################################
# Update Stress Last
######################################################################

function solve_explicit_dynamics_femp_3D(grid,solids,basis,body,alg::USL,output,fixes,data)
   
  Tf             = data["total_time"]::Float64
	dtime          = data["dt"]        ::Float64     
	t              = data["time"]      ::Float64
  counter        = 0 

  Identity       = UniformScaling(1.)
	solidCount     = length(solids)
	nodalMass      = grid.mass
	nodalMomentum0 = grid.momentum0
	nodalMomentum  = grid.momentum
	nodalForce     = grid.force

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
    D         = SMatrix{3,3}(0,0,0,0,0,0,0,0,0) #zeros(Float64,2,2)
    g         = [0.,0.,0.]

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
			#du[ip] = @SVector [0., 0., 0.]
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
	
		   	D           = 0.5 * (vel_grad + vel_grad')
		   	strain[ip]   = D
		   	F[ip]        = inv(Identity - vel_grad)
		   	J            = det(F[ip])
		   	#println(strain[ip])
	   	     #@timeit "3" update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
	   	    update_stress!(stress[ip],mat,strain[ip],D,F[ip],J,ip,dtime)

            sigma = stress[ip]
            body(g,xx[ip],t)  
            #println(g)
            # compute nodal internal force fint
		    for i = 1:length(elemNodes)
				  in  = elemNodes[i]; # index of node 'i'
			    dNi = @view dNdx[:,i]			
	   	        fint[in]+= w*@SVector[sigma[1,1] * dNi[1] + sigma[1,2] * dNi[2] + sigma[1,3] * dNi[3],
									  sigma[1,2] * dNi[1] + sigma[2,2] * dNi[2] + sigma[2,3] * dNi[3],
									  sigma[1,3] * dNi[1] + sigma[2,3] * dNi[2] + sigma[3,3] * dNi[3] ]
                fbody[in] += w*N[i]*mat.density*g										  
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

# data is a dictionary: 
# data["pressure"]   = [(1,"tag1",f),(2,"tag2",g)]
# data["total_time"] = Tf
# data["dt"]         = dtime
# data["time"]       = t


function solve_explicit_dynamics_femp_3D(grid,solids,mats,basis,body,alg::TLFEM,output,fixes,data)
  Tf    = data["total_time"]::Float64
	dtime = data["dt"]        ::Float64     
	t     = data["time"]      ::Float64

  counter = 0

  Identity       = UniformScaling(1.)
	solidCount     = length(solids)
	nodalMass      = grid.mass
	nodalMomentum0 = grid.momentum0
	nodalMomentum  = grid.momentum
	nodalForce     = grid.force

    # pre_allocating arrays for temporary variable
   
	nearPoints,funcs, ders = initialise(grid,basis)

  nodePerElem = size(solids[1].elems,2)

	if nodePerElem == 4 		
		# wgt       = .166666667
		# weights   = [0.166666667]   # this is so dangerous!!! all books ay weight =1
		# gpCoords  = [0.25,0.25,0.25]
		
        weights_surface   = ones(3)
        normals_surface   = zeros(3, 1)
        funcs_surface     = zeros(3,1)
        gpCoords_surface  = zeros(2,1)
	end
	if nodePerElem == 8 		
		# wgt       = 8.0
  #       weights   = [8.]
  #       gpCoords  = [0.,0.,0.]

        # for pressure load
        weights_surface   = ones(4)
        normals_surface   = zeros(3, 4)
        funcs_surface     = zeros(4,4)
        gpCoords_surface  = zeros(2,4)

        gpCoords_surface[1,1] = -0.5773502691896257; gpCoords_surface[2,1] = -0.5773502691896257;
        gpCoords_surface[1,2] =  0.5773502691896257; gpCoords_surface[2,2] = -0.5773502691896257;
        gpCoords_surface[1,3] =  0.5773502691896257; gpCoords_surface[2,3] =  0.5773502691896257;
        gpCoords_surface[1,4] = -0.5773502691896257; gpCoords_surface[2,4] =  0.5773502691896257;
	end 

    vel_grad  = SMatrix{3,3}(0,0,0,0,0,0,0,0,0)
    D         = SMatrix{3,3}(0,0,0,0,0,0,0,0,0) #zeros(Float64,2,2)
    g         = [0.,0.,0.]
    Fdot      = SMatrix{3,3}(0,0,0,0,0,0,0,0,0)
   

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
        initializeBasis(solid,mats[s].density)
    end

    # time-independent Dirichlet boundary conditions on grid/solids
    fix_Dirichlet_grid(grid,data)
    fix_Dirichlet_solid(solids,data)

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
	
	    # deformable solids only
		if typeof(mats[s]) <: RigidMaterial continue end

		xx     = solid.pos
		mm     = solid.mass
		vv     = solid.velocity		
		fint   = solid.fint
		fbody  = solid.fbody
		ftrac  = solid.ftrac
		du     = solid.dU

	  	@inbounds for ip = 1:solid.nodeCount	        
			support   = getShapeFunctions(nearPoints,funcs,ip, grid, solid,basis)
	        
	        fMass     = mm[ip]
	        vp        = vv[ip]
	        fp        = fint[ip]
	        fb        = fbody[ip]
	        ft        = ftrac[ip]
			@inbounds for i = 1:support
				in    = nearPoints[i]; # index of node 'i'
				Ni    = funcs[i]				
				Nim   = Ni * fMass
				# mass, momentum, internal force and external force
		    nodalMass[in]      += Nim
				nodalMomentum0[in] += Nim * vp 
		    nodalForce[in]     -= Ni  * fp
				nodalForce[in]     += Ni  * fb
				nodalForce[in]     += Ni  * ft
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
        fixed_dirs       = @view grid.fixedNodes[:,i]
        if fixed_dirs[1] == 1
			nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,1)
   		    nodalMomentum[i]  = setindex(nodalMomentum[i], 0.,1)
        end
        if fixed_dirs[2] == 1
			nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,2)
			nodalMomentum[i]  = setindex(nodalMomentum[i], 0.,2)
        end
        if fixed_dirs[3] == 1
			nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,3)
			nodalMomentum[i]  = setindex(nodalMomentum[i], 0.,3)
        end
	end

	# ===========================================
	# particle to grid (rigid solids)
	# ===========================================

	if haskey(data, "rigid_body_velo") 

	    rigid_solids = data["rigid_body_velo"]::Array{Tuple{Int64,Function},1}

	 	for (s,f) in rigid_solids
	 		solid       = solids[s]
			xx          = solid.pos
		
			vex,vey,vez = f(t)::Tuple{Float64,Float64,Float64}
			
			@inbounds for ip = 1:solid.nodeCount
				getAdjacentGridPoints(nearPoints,xx[ip],grid,basis)
				#support   = getShapeAndGradient(nearPoints,funcs,ip, grid, solid,basis)
	#			println(nearPoints)
				@inbounds for i = 1:8
					id                  = nearPoints[i]; # index of node 'i'
					mi                  = nodalMass[id]
	
				    nodalMomentum0[id]  = mi *@SVector[vex,vey,vez]	
				    nodalMomentum[id]   = mi *@SVector[vex,vey,vez]	

				end
			end
			#solid.reaction_forces .= @SVector[Fx, Fy, Fz]
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

			#du[ip] = setindex(du[ip],0.,3)
			#vv[ip] = setindex(vv[ip],0.,3)

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

    fix_Dirichlet_solid(solids,data,dtime)

    # ====================================================================
    #  update particle internal forces
    # ====================================================================

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
	  	neighbours = solid.neighbours
	  	
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
			
			Fdot         = vel_grad/dtime

		   	
		   	F[ip]        = Identity + vel_gradT
		   	J            = det(F[ip])

		   	vol[ip]      = detJ*J

		   	L             = Fdot*inv(F[ip])
		   	D             = 0.5 * dtime * (L + L')
		   	strain[ip]   += D #0.5 * (vel_grad + vel_grad' + vel_grad * vel_grad')

		   	#println(strain[ip])
	   	     #@timeit "3" update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)

	   	    stress[ip] =  update_stress!(stress[ip],mat,strain[ip],D,F[ip],J,ip,dtime)


            P     = J*stress[ip]*inv(F[ip])'  # convert to Piola Kirchoof stress

            body(g,xx[ip],t)  
            #println(g)
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

	# ====================================================================
    #  update particle external forces due to traction
    # ====================================================================

    #t       += dtime

    compute_fext(solids,funcs_surface, normals_surface, weights_surface, gpCoords_surface,data,t)

    # checking the pressure 
 #    ftrac = solids[1].ftrac
 #    sss = 0.
	# for ip=1:solids[1].nodeCount 
 #      sss += sqrt(ftrac[ip][1]^2+ftrac[ip][2]^2+ftrac[ip][3]^2)
	# end
	# println(sss)

    if haskey(data, "rigid_body_velo") 
	    
	    rigid_solids = data["rigid_body_velo"]::Array{Tuple{Int64,Function},1}

	 		for (s,f) in rigid_solids
		 		solid       = solids[s]
				xx          = solid.pos
				#xc          = solid.centroids
				(vex,vey,vez) = f(t)::Tuple{Float64, Float64, Float64} 

				if (vex,vey,vez) == (0.,0.,0.) continue end
				# update all nodes for visualization only
				@inbounds for ip = 1:solid.nodeCount
			      xx[ip]   += dtime * @SVector [vex,vey,vez]   # no memory alloc, 
			  end
	            # update the centroids of surface elems
			    #@inbounds for ip = 1:solid.surfCount
			    #  xc[ip]   += dtime * @SVector [vex,vey,vez]
			    #end
		  end
	  end

    # write stuff to files at intervals
		if (counter%output.interval == 0)
			#println("haha\n")
			plotParticles_3D(output,solids,mats,counter)
			#plotGrid(output,grid,counter)
			compute_femp(fixes,t)
		end

    # advance to next time step
    t       += dtime
    counter += 1
  end # end of time loop
  #@printf("Solving done \n")
  closeFile(fixes)
end # end solve()


function solve_explicit_dynamics_femp_eigen_erosion_3D(grid,solids,mats,basis,body,alg::TLFEM,output,fixes,data)
  Tf    = data["total_time"]
	dtime = data["dt"]         
	t     = data["time"]      

  counter = 0

  Identity       = UniformScaling(1.)
	solidCount     = length(solids)
	nodalMass      = grid.mass
	nodalMomentum0 = grid.momentum0
	nodalMomentum  = grid.momentum
	nodalForce     = grid.force

    # pre_allocating arrays for temporary variable
   
	nearPoints,funcs, ders = initialise(grid,basis)

    nodePerElem = size(solids[1].elems,2)

	if nodePerElem == 4 		
		# wgt       = .166666667
		# weights   = [0.166666667]   # this is so dangerous!!! all books ay weight =1
		# gpCoords  = [0.25,0.25,0.25]


        # for pressure load
        weights_surface   = ones(3)
        normals_surface   = zeros(3, 1)
        funcs_surface     = zeros(3,1)
        gpCoords_surface  = zeros(2,1)

	end
	if nodePerElem == 8 		
		# wgt       = 8.0
  #       weights   = [8.]
  #       gpCoords  = [0.,0.,0.]

        # for pressure load
        weights_surface   = ones(4)
        normals_surface   = zeros(3, 4)
        funcs_surface     = zeros(4,4)
        gpCoords_surface  = zeros(2,4)

        gpCoords_surface[1,1] = -0.5773502691896257; gpCoords_surface[2,1] = -0.5773502691896257;
        gpCoords_surface[1,2] =  0.5773502691896257; gpCoords_surface[2,2] = -0.5773502691896257;
        gpCoords_surface[1,3] =  0.5773502691896257; gpCoords_surface[2,3] =  0.5773502691896257;
        gpCoords_surface[1,4] = -0.5773502691896257; gpCoords_surface[2,4] =  0.5773502691896257;
	end 

    vel_grad  = SMatrix{3,3}(0,0,0,0,0,0,0,0,0)
    D         = SMatrix{3,3}(0,0,0,0,0,0,0,0,0) #zeros(Float64,2,2)
    g         = [0.,0.,0.]
    Fdot      = zeros(3,3)
   

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
        initializeBasis(solid,mats[s].density)
    end

    # time-independent Dirichlet boundary conditions on grid/solids
    fix_Dirichlet_grid(grid,data)
    fix_Dirichlet_solid(solids,data)

    α = data["alpha"]
    κ = data["kappa"]
    build_particle_neighbors(solids,grid,κ)

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
		ftrac  = solid.ftrac
		du     = solid.dU

	  	@inbounds for ip = 1:solid.nodeCount	        
			support   = getShapeFunctions(nearPoints,funcs,ip, grid, solid,basis)
	        
	        fMass     = mm[ip]
	        vp        = vv[ip]
	        fp        = fint[ip]
	        fb        = fbody[ip]
	        ft        = ftrac[ip]
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
				nodalForce[in]     += Ni  * ft
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
        fixed_dirs       = @view grid.fixedNodes[:,i]
        if fixed_dirs[1] == 1
			nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,1)
   		    nodalMomentum[i]  = setindex(nodalMomentum[i], 0.,1)
        end
        if fixed_dirs[2] == 1
			nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,2)
			nodalMomentum[i]  = setindex(nodalMomentum[i], 0.,2)
        end
        if fixed_dirs[3] == 1
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
	  	fix   = solid.fixedNodes
	  	
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

			#du[ip] = setindex(du[ip],0.,3)
			#vv[ip] = setindex(vv[ip],0.,3)

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

    fix_Dirichlet_solid(solids,data,dtime)

    # ====================================================================
    #  update particle internal forces
    # ====================================================================

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
	  	neighbours = solid.neighbours
	  	elem_mass = solid.elem_mass 
	  	
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
		   	F[ip]        = Identity + vel_gradT
		   	L             = Fdot*inv(F[ip])
		   	D             = 0.5 * dtime * (L + L')
		   	strain[ip]   += D #0.5 * (vel_grad + vel_grad' + vel_grad * vel_grad')

		   	#println(strain[ip])
	   	     #@timeit "3" update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)

	   	    update_stress!(stress[ip],mat,strain[ip],D,F[ip],det(F[ip]),ip,dtime)
	   end    
       
        # ################################################
        # # eigen-erosion farcture
        # ################################################
	  	@inbounds for ip = 1:solid.parCount
		
            totalH = 0.
            totalM = 0.

            for c in neighbours[ip]
            	# skip eroded particle 
            	if mat.dam[c] >= 1. continue end 
            	H  = computeCrackDrivingForce(stress[c],mat)
            	mc = elem_mass[c]
            	totalH += mc * H
            	totalM += mc 
            end

            energy = κ * α * (totalH / totalM)
            if (energy >= mat.Gf) 
            	mat.dam[ip] = 1.0 
            else
            	J     = det(F[ip])
	            P     = J*stress[ip]*inv(F[ip])'  # convert to Piola Kirchoof stress

	            #body(g,[0 0 0],t)  
	            elemNodes =  @view elems[ip,:]  
	            dNdx      = feDers[ip]
	            # compute nodal internal force fint
			    for i = 1:length(elemNodes)
					in  = elemNodes[i]; # index of node 'i'
				    dNi = @view dNdx[:,i]			
		   	        fint[in]  +=  vol[ip] * @SVector[P[1,1] * dNi[1] + P[1,2] * dNi[2] + P[1,3] * dNi[3],
											      P[2,1] * dNi[1] + P[2,2] * dNi[2] + P[2,3] * dNi[3],
											      P[3,1] * dNi[1] + P[3,2] * dNi[2] + P[3,3] * dNi[3] ]	
	                #fbody[in] += vol[ip]  * mat.density*N[i]*g								        
	            end
            end
	   end

	end

	# ====================================================================
    #  update particle external forces due to traction
    # ====================================================================

    t       += dtime

    compute_fext(solids,funcs_surface, normals_surface, weights_surface, gpCoords_surface,data,t)

    # checking the pressure 
 #    ftrac = solids[1].ftrac
 #    sss = 0.
	# for ip=1:solids[1].nodeCount 
 #      sss += sqrt(ftrac[ip][1]^2+ftrac[ip][2]^2+ftrac[ip][3]^2)
	# end
	# println(sss)

	if (counter%output.interval == 0)
		#println("haha\n")
		plotParticles_3D(output,solids,mats,counter)
		#plotGrid(output,grid,counter)
		compute_femp(fixes,t)
	end

    
    counter += 1
  end # end of time loop
  #@printf("Solving done \n")
  closeFile(fixes)
end # end solve()


######################################################################
# Update Stress Last, TLFEM for internal force, axi-symmetric
######################################################################
function solve_explicit_dynamics_femp_3D(grid,solids::Array{FEMAxis,1},mats,basis,body,alg::TLFEM,output,fixes,Tf,dtime)
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

    vel_grad   = zeros(3,3)
    vel_gradT  = zeros(3,3)
    D          = SMatrix{3,3}(0,0,0,0,0,0,0,0,0) #zeros(Float64,2,2)
    g          = [0.,0.]
  

	@inbounds for s = 1:solidCount
		solid  = solids[s]
        initializeBasis(solid,mats[s].density)
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
		ftrac  = solid.ftrac
		du     = solid.dU

	  	@inbounds for ip = 1:solid.nodeCount	        
			support   = getShapeFunctions(nearPoints,funcs,ip, grid, solid,basis)
	        
	        fMass     = mm[ip]
	        vp        = vv[ip]
	        fp        = fint[ip]
	        fb        = fbody[ip]
	        ft        = ftrac[ip]
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
				nodalForce[in]     += Ni  * ft
			end
			du[ip] = @SVector [0., 0.]
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
	  		fint[ip]  = @SVector [0., 0.]
	  		fbody[ip] = @SVector [0., 0.]
	  	end
        # loop over solid elements, not solid nodes
	  	@inbounds for ip = 1:solid.parCount
			elemNodes =  @view elems[ip,:]  			
			detJ      = feJaco[ip]
			N         = feFuncs[ip]
			dNdx      = feDers[ip]
			
			vol[ip]   = detJ
			#println(dNdx)
			#println(sum(dNdx, dims=2))
			vel_grad  .= 0.
			vel_gradT .= 0.
			for i = 1:length(elemNodes)
				in         = elemNodes[i]; # index of node 'i'
			    dNi        = @view dNdx[:,i]			
				vI         = du[in]
				vIt        = xx[in] - XX[in]

				vel_grad[1,1]  += dNi[1]*vI[1] 
				vel_grad[1,2]  += dNi[2]*vI[1] 
				vel_grad[2,1]  += dNi[1]*vI[2] 
				vel_grad[2,2]  += dNi[2]*vI[2] 
				vel_grad[3,3]  += N[i]*vI[1] 

				vel_gradT[1,1]  += dNi[1]*vIt[1] 
				vel_gradT[1,2]  += dNi[2]*vIt[1] 
				vel_gradT[2,1]  += dNi[1]*vIt[2] 
				vel_gradT[2,2]  += dNi[2]*vIt[2] 
				vel_gradT[3,3]  += N[i]*vIt[1] 
		   	end
		   	
			#dstrain      = 0.5 * (vel_grad + vel_grad' + vel_grad * vel_grad') - strain[ip] 
			
			Fdot         = vel_grad/dtime

		   	
		   	F[ip]        = Identity + vel_gradT
		   	J            = det(F[ip])

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
	   	        fint[in]  +=  detJ * @SVector[P[1,1] * dNi[1] + P[1,2] * dNi[2] + P[3,3] * N[i],
										      P[2,1] * dNi[1] + P[2,2] * dNi[2]  ]
                fbody[in] += detJ*mat.density*N[i]*g								        
            end
	   end
	end

	t       += dtime

	if (counter%output.interval == 0)
		#println("haha\n")
		plotParticles_2D(output,solids,mats,counter)
		#plotGrid(output,grid,counter)
		compute_femp(fixes,t)
	end

    
    counter += 1
  end # end of time loop
  #@printf("Solving done \n")
  closeFile(fixes)
end # end solve()


######################################################################
# Update Stress Last, TLFEM  with full quadrature for internal force
######################################################################

function solve_explicit_dynamics_femp_3D(grid,solids,mats, basis,body,alg::TLFEMFull,output,fixes,data)
  Tf             = data["total_time"]::Float64
	dtime          = data["dt"]        ::Float64     
	t              = data["time"]      ::Float64
  counter        = 0

  Identity       = UniformScaling(1.)
	solidCount     = length(solids)
	nodalMass      = grid.mass
	nodalMomentum0 = grid.momentum0
	nodalMomentum  = grid.momentum
	nodalForce     = grid.force

  D              = SMatrix{3,3}(0,0,0,0,0,0,0,0,0) #zeros(Float64,2,2)

    # pre_allocating arrays for temporary variable
   
	nearPoints,funcs, ders = initialise(grid,basis)

  nodePerElem = size(solids[1].elems,2)

  # four-node Tetrahedron elements
	if nodePerElem == 4 
		meshBasis = Tet4()  
		wgt       = 1.
		weights   = [1.]
		gpCoords  = [0.25,0.25,0.25]


    # for pressure load
    weights_surface   = ones(4)
    normals_surface   = zeros(3, 4)
    funcs_surface     = zeros(4,4)
	end

	# eight-node hexahedron elements
	if nodePerElem == 8 
		meshBasis = Hexa8() 
		#wgt       = 8.0
    #weights   = [8.]
    #gpCoords  = [0.,0.,0.]

	  gpCoords  = zeros(3,8)
	  weights   = ones(8)
	  gpCoords[1,1] = -0.5773502691896257;gpCoords[2,1] = -0.5773502691896257;gpCoords[3,1] = -0.5773502691896257
	  gpCoords[1,2] =  0.5773502691896257;gpCoords[2,2] = -0.5773502691896257;gpCoords[3,2] = -0.5773502691896257
	  gpCoords[1,3] =  0.5773502691896257;gpCoords[2,3] =  0.5773502691896257;gpCoords[3,3] = -0.5773502691896257
	  gpCoords[1,4] = -0.5773502691896257;gpCoords[2,4] =  0.5773502691896257;gpCoords[3,4] = -0.5773502691896257
	  gpCoords[1,5] = -0.5773502691896257;gpCoords[2,5] =  0.5773502691896257;gpCoords[3,5] =  0.5773502691896257
	  gpCoords[1,6] = -0.5773502691896257;gpCoords[2,6] =  0.5773502691896257;gpCoords[3,6] =  0.5773502691896257
	  gpCoords[1,7] = -0.5773502691896257;gpCoords[2,7] =  0.5773502691896257;gpCoords[3,7] =  0.5773502691896257
	  gpCoords[1,8] = -0.5773502691896257;gpCoords[2,8] =  0.5773502691896257;gpCoords[3,8] =  0.5773502691896257

    # for pressure load
    weights_surface   = ones(4)
    normals_surface   = zeros(3, 4)
    funcs_surface     = zeros(4,4)
    gpCoords_surface  = zeros(2,4)

    gpCoords_surface[1,1] = -0.5773502691896257; gpCoords_surface[2,1] = -0.5773502691896257;
    gpCoords_surface[1,2] =  0.5773502691896257; gpCoords_surface[2,2] = -0.5773502691896257;
    gpCoords_surface[1,3] =  0.5773502691896257; gpCoords_surface[2,3] =  0.5773502691896257;
    gpCoords_surface[1,4] = -0.5773502691896257; gpCoords_surface[2,4] =  0.5773502691896257;
	end

  noGP      = size(gpCoords,2)
  dNdx      = zeros(3,nodePerElem)
  N         = zeros(nodePerElem)#@SVector [0,0,0,0]

  vel_grad  = SMatrix{3,3}(0,0,0,0,0,0,0,0,0)
  g         = [0.,0.,0.]
   

  # compute nodal mass (only once)


  @inbounds for s = 1:solidCount
		solid  = solids[s]
		xx     = solid.pos
		mm     = solid.mass  # to be updated here
		elems  = solid.elems
		mat    = mats[s]
		rho    = mats[s].density
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
			  	mm[id]  += rho*N[i]*detJ*weights[gp]
			  	vol[ip] += detJ*weights[gp]
			  end
			end
	  	end
    end

  # time-independent Dirichlet boundary conditions on grid/solids
  fix_Dirichlet_grid(grid,data)
  fix_Dirichlet_solid(solids,data)

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
		ftrac  = solid.ftrac
		du     = solid.dU

	  	@inbounds for ip = 1:solid.nodeCount	        
			support   = getShapeFunctions(nearPoints,funcs,ip, grid, solid,basis)
	        
	        fMass     = mm[ip]
	        vp        = vv[ip]
	        fp        = fint[ip]
	        fb        = fbody[ip]
	        ft        = ftrac[ip]
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
				nodalForce[in]     += Ni  * ft
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
    fixed_dirs       = @view grid.fixedNodes[:,i]
    if fixed_dirs[1] == 1
			nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,1)
   		    nodalMomentum[i]  = setindex(nodalMomentum[i], 0.,1)
    end
    if fixed_dirs[2] == 1
			nodalMomentum0[i] = setindex(nodalMomentum0[i],0.,2)
			nodalMomentum[i]  = setindex(nodalMomentum[i], 0.,2)
    end
    if fixed_dirs[3] == 1
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

      if typeof(mats[s]) <: RigidMaterial continue end

	  	xx    = solid.pos
	  	mm    = solid.mass
	  	vv    = solid.velocity
	  	du    = solid.dU
	  	fix   = solid.fixedNodes
	  	
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

			#du[ip] = setindex(du[ip],0.,3)
			#vv[ip] = setindex(vv[ip],0.,3)

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

    fix_Dirichlet_solid(solids,data,dtime)
	


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
	  	mat    = mats[s]
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
	   	    stress[ip] =  update_stress!(stress[ip],mat,strain[ip],D,F[ip],J,ip,dtime)

            sigma = stress[ip]
            P     = J*sigma*inv(F[ip])'  # convert to Piola Kirchoof stress

            body(g,XX[ip],t)  
            # compute nodal internal force fint
		    for i = 1:length(elemNodes)
				in  = elemNodes[i]; # index of node 'i'
			    dNi = @view dNdx[:,i]			
	   	        fint[in]  +=  w * @SVector[P[1,1] * dNi[1] + P[1,2] * dNi[2] + P[1,3] * dNi[3],
										   P[2,1] * dNi[1] + P[2,2] * dNi[2] + P[2,3] * dNi[3],
										   P[3,1] * dNi[1] + P[3,2] * dNi[2] + P[3,3] * dNi[3] ]
                fbody[in] += w*mat.density*N[i]*g								        
            end
	   end
	end

	# ====================================================================
    #  update particle external forces due to traction
    # ====================================================================

    t       += dtime

    compute_fext(solids,funcs_surface, normals_surface, weights_surface, gpCoords_surface,data,t)

 #    @inbounds for s = 1:solidCount
	# 	# only deformable solids here
	#   	solid = solids[s]

	#   	XX     = solid.pos0	  	  	  
	#   	elems  = solid.mesh.elements
	#   	ftrac  = solid.ftrac
	#   	for ip=1:solid.nodeCount 
	#   		ftrac[ip]  = @SVector [0., 0., 0.]	  		
	#   	end

 #        if haskey(solid.mesh.element_sets, "force")
	#   	    surf_elems_ids = collect(solid.mesh.element_sets["force"])
	#   	else
	#   		continue
	#   	end
 #        # loop over  elements of the surface tag 'force'
	#   	@inbounds for ip in surf_elems_ids
	# 		elemNodes =  elems[ip]  
	# 		coords    =  @view XX[elemNodes]
	# 		#println(coords)
	# 		getNormals!(funcs_surface, normals_surface, weights_surface , coords, gpCoords_surface, Quad4() )
	# 		# loop over Gauss points
	# 		for ip=1:4                     
	# 		    ww = weights_surface[ip]   
	# 		    #println(normals_surface[:,ip])       
	# 		    #println(ww)       
	# 			for i = 1:length(elemNodes)
	# 			   in  = elemNodes[i]; # index of node 'i'
	# 	           ftrac[in][1] -= normals_surface[1,ip]*funcs_surface[i,ip]*ww*400*exp(-10000*t)
	# 	           ftrac[in][2] -= normals_surface[2,ip]*funcs_surface[i,ip]*ww*400*exp(-10000*t)
	# 	           ftrac[in][3] -= normals_surface[3,ip]*funcs_surface[i,ip]*ww*400*exp(-10000*t)
	# 	        end
	#         end
	#    end
	# end


	if (counter%output.interval == 0)
		#println("haha\n")
		plotParticles_3D(output,solids,mats,counter)
		#plotGrid(output,grid,counter)
		compute_femp(fixes,t)
	end

    
    counter += 1
  end # end of time loop
  #@printf("Solving done \n")
  closeFile(fixes)
end # end solve()
