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
push!(LOAD_PATH,"./")
using Velocity
using Printf
######################################################################
# Update Stress Last, TLFEM for internal force
######################################################################
function solve_explicit_dynamics_femp_3D_Contact(grid,solids,mats,basis,body,alg::TLFEM,output,fixes,data)
    Tf    = data["total_time"]    :: Float64
    dtime = data["dt"]            :: Float64
    t     = data["time"]          :: Float64
    fric  = data["friction"]      :: Float64
    dt_factor = data["dt_factor"] :: Float64
    counter = 0

    Identity       = UniformScaling(1.)
    solidCount     = length(solids)

    # this is body mass/momenta/forces/normals
    nodalMass      =  [copy(grid.mass)       for _ in 1:solidCount]
    nodalMomentum0 =  [copy(grid.momentum0)  for _ in 1:solidCount]
    nodalMomentum  =  [copy(grid.momentum)   for _ in 1:solidCount]
    nodalForce     =  [copy(grid.force)      for _ in 1:solidCount]
    nodalNormals   =  [copy(grid.force)      for _ in 1:solidCount]
    nodalPos       =  [copy(grid.force)      for _ in 1:solidCount]  # used for better contact checking [Nairn]

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
    contact_solids   = Vector{Vector{Int64}}(undef,nodeCount)

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
    deltaVe1  = @SVector [0., 0., 0.] #zeros(3)
    dVxNI     = @SVector [0., 0., 0.] #zeros(3)
    omega     = @SVector [0., 0., 0.] #zeros(3)
    nI        = @SVector [0., 0., 0.] #zeros(3)
    nI_common = @SVector [0., 0., 0.] #zeros(3)

    # compute nodal mass and FE shape functions and derivatives at centroids for all
    # elements in the mesh of each solid
    @inbounds for s = 1:solidCount
	solid  = solids[s]
	#if typeof(mats[s]) <: RigidMaterial
	compute_normals!(solid)
	#	continue
	#end
        initializeBasis(solid,mats[s].density)
    end

    # time-independent Dirichlet boundary conditions on grid/solids
    fix_Dirichlet_grid(grid,data)
    fix_Dirichlet_solid(solids,data)

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


    #################################################
    # Time loop 
    #################################################

    while t < Tf
        #if counter > 10
        #    break
        #end

        @printf("Solving step: %d %e %f\n", counter, dtime, t)

        # ===========================================
        # reset grid data
        # ===========================================

        @inbounds for i = 1:nodeCount
	    nodalMass_S[i]      = 0.	  
	    nodalMomentum_S[i]  =  @SVector [0., 0., 0.]
            if counter == 0
                contact_solids[i] = Vector{Int64}()#Set{Int64}()
            else
                #@inbounds for s = 1:solidCount
                #    delete!(contact_solids[i], s)
                #end
                empty!(contact_solids[i])
            end
        end

        boundary_nodes = Set{Int64}()

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

	    nMass     = nodalMass[s]
	    nMomenta0 = nodalMomentum0[s]		
	    nMomenta  = nodalMomentum[s]		
	    nForce    = nodalForce[s]
	    normals   = nodalNormals[s]
	    x_node    = nodalPos[s]
	    bnd_particles = solid.mesh.node_sets["boundary"]	

            # reset body mass/momenta0/force 
	    nMass[:]  .= 0.
	    fill!(nMomenta0,@SVector[0.,0.,0.])
	    fill!(nForce,   @SVector[0.,0.,0.])
	    fill!(normals,  @SVector[0.,0.,0.])
	    fill!(x_node,   @SVector[0.,0.,0.])
	    #bnd_nodes     = Set{Int64}() # old implementation: solid.boundary_nodes
	    
	    @inbounds for ip = 1:solid.nodeCount
	        #getShapeAndGradient(nearPoints,funcs,ders,xx[ip], grid) 
		support   = getShapeAndGradient(nearPoints,funcs,ders,ip, grid, solid,basis)
	        fMass     = mm[ip]
	        vp        = vv[ip]
	        fp        = fint[ip]
	        fb        = fbody[ip]
	        ft        = ftrac[ip]
		#println(funcs)
		@inbounds for i = 1:support
		    id    = nearPoints[i]; # index of node 'i'
		    Ni    = funcs[i]				
		    dNi   = @view ders[:,i]			
		    Nim   = Ni * fMass
		    # mass, momentum, internal force and external force
                    nMass[id]      += Nim             # check: no memory alloc
                    nMomenta0[id]  += Nim * vp
		    nForce[id]     -= Ni  * fp
                    nForce[id]     += Ni  * fb       # check: no memory alloc
                    nForce[id]     += Ni  * ft
		end
		du[ip] = @SVector [0., 0., 0.]
	    end

            ############################################################################
	    # another loop over surface nodes => grid boundary nodes and normals
            ######################################################################
	    #xx          = solid.centroids		
	    surf_normals = solid.normals

	    # loop over centroids of surface elements
	    @inbounds for ip in collect(solid.mesh.node_sets["boundary"])
		support    = getShapeFunctions(nearPoints,funcs,xx[ip],grid, basis)
		#println(nearPoints)
		normal_ip = surf_normals[ip]
		xp        = xx[ip]
		mp        = mm[ip]
                #			println(nearPoints)
		@inbounds for i = 1:8
		    id                  = nearPoints[i]; # index of node 'i'
		    Ni                  = funcs[i]					
                    normals[id]        += Ni * normal_ip
                    x_node[id]         += Ni * xp * mp
		    #if (ip in bnd_particles ) push!(boundary_nodes,id) end
                    if !(s in contact_solids[id])
		        push!(contact_solids[id],s)  # add solid 's' to the set of contact solids at node 'id'
                    end
		    push!(boundary_nodes,id) 
		end
	    end # end loop over centroids


	    #boundary_nodes[s] = bnd_nodes old implementation
	    
	    # ===========================================
	    # update grid
	    # ===========================================

	    @inbounds for i=1:nodeCount

                nMomenta[i]         = nMomenta0[i] + dtime * nForce[i]
                nodalMass_S[i]      += nMass[i]
		nodalMomentum_S[i]  += nMomenta[i]

                nMomenta0[i]  /= nMass[i]
		nMomenta[i]   /= nMass[i]
                #if i == 47488
                #    @printf("nMomenta[%d] = [%f, %f, %f]\n", i, nMomenta[i][1], nMomenta[i][2], nMomenta[i][3])
                #end
		

	        # apply Dirichet boundary conditions on which solid?

	        
	        fixed_dirs       = @view grid.fixedNodes[:,i]
	        if fixed_dirs[1] == 1
                    nMomenta[i]  = setindex(nMomenta[i],0.,1) 
                    nMomenta0[i] = setindex(nMomenta0[i],0.,1)# => memory alloc!!!	    
	        end
	        if fixed_dirs[2] == 1
		    nMomenta[i]  = setindex(nMomenta[i],0.,2) 
                    nMomenta0[i] = setindex(nMomenta0[i],0.,2)# => memory alloc!!!
	        end
	        if fixed_dirs[3] == 1
		    nMomenta[i]  = setindex(nMomenta[i],0.,3) 
                    nMomenta0[i] = setindex(nMomenta0[i],0.,3)# => memory alloc!!!
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
		xx          = solid.pos		
		surf_normals= solid.normals
		normals     = nodalNormals[s]
		x_node      = nodalPos[s]
	        
		fill!(normals,  @SVector[0.,0.,0.])
		fill!(x_node,   @SVector[0.,0.,0.])
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
                        normals[id]        += funcs[i] * normal_ip
			mi                  = nodalMass_S[id]
			#if (ip in bnd_particles ) push!(boundary_nodes,id) end
                        if !(s in contact_solids[id])
			    push!(contact_solids[id],s)  # add solid 's' to the set of contact solids at node 'id'
                        end

                        # Fx += (mi*vex - nodalMomentum_S[id][1])/dtime
                        # Fy += (mi*vey - nodalMomentum_S[id][2])/dtime
                        # Fz += (mi*vez - nodalMomentum_S[id][3])/dtime
	                
                        nodalMomentum_S[id]   =  mi*@SVector[vex,vey,vez]
                        #x_node[id]            = grid.pos[id] + t * @SVector[vex,vey,vez]
                        #if id==47488
                        #    dN = funcs[i] * normal_ip
                        #    @printf("solid.surfCount=%d\n",solid.surfCount)
                        #    @printf("xx[%d] = [%f, %f, %f], funcs[i]=%.3e\n", ip, xx[ip][1], xx[ip][2], xx[ip][3], funcs[i])
                        #    @printf("normal_ip = [%.3e, %.3e, %.3e]\n", normal_ip[1], normal_ip[2], normal_ip[3])
                        #    @printf("dN = [%.3e, %.3e, %.3e]\n", dN[1], dN[2], dN[3])
                        #    @printf("normals[%d] = [%.3e, %.3e, %.3e]\n", id, normals[id][1], normals[id][2], normals[id][3])
                        #end
		    end
		end
		#solid.reaction_forces .= @SVector[Fx, Fy, Fz]
	    end
	end

	if haskey(data, "rigid_body_velo_data") 

	    rigid_solids = data["rigid_body_velo_data"]::Array{Tuple{Int64,VelocityData},1}

	    for (s,f) in rigid_solids
	 	solid       = solids[s]
		xx          = solid.pos
		surf_normals= solid.normals
		normals     = nodalNormals[s]
		x_node      = nodalPos[s]
	        
		fill!(normals,  @SVector[0.,0.,0.])
		fill!(x_node,   @SVector[0.,0.,0.])
		(vex,vey,vez) = f(t)::Tuple{Float64, Float64, Float64} 
                @printf("\tve = [%4.3e, %4.3e, %4.3e]\t", vex, vey, vez)
		#Fx = Fy = Fz = 0.
		# loop over centroids of surface elements
		@inbounds for ip = 1:solid.surfCount
		    support    = getShapeFunctions(nearPoints,funcs,xx[ip],grid, basis)
		    #println(nearPoints)
		    normal_ip = surf_normals[ip]
	            #			println(nearPoints)
		    @inbounds for i = 1:8
			id                  = nearPoints[i]; # index of node 'i'					
                        normals[id]        += funcs[i] * normal_ip
			mi                  = nodalMass_S[id]
			#if (ip in bnd_particles ) push!(boundary_nodes,id) end
                        if !(s in contact_solids[id])
			    push!(contact_solids[id],s)  # add solid 's' to the set of contact solids at node 'id'
                        end

                        # Fx += (mi*vex - nodalMomentum_S[id][1])/dtime
                        # Fy += (mi*vey - nodalMomentum_S[id][2])/dtime
                        # Fz += (mi*vez - nodalMomentum_S[id][3])/dtime
	                
                        nodalMomentum_S[id]   =  mi*@SVector[vex,vey,vez]
                        #if mi != 0
                        #    @printf("nodalMomentum_S[%d] = [%4.3e, %4.3e, %4.3e]\n", id, nodalMomentum_S[id][1], nodalMomentum_S[id][2], nodalMomentum_S[id][3])
                        #end
                        #x_node[id]            = grid.pos[id] + t * @SVector[vex,vey,vez]
                        #if id==47488
                        #    dN = funcs[i] * normal_ip
                        #    @printf("solid.surfCount=%d\n",solid.surfCount)
                        #    @printf("xx[%d] = [%f, %f, %f], funcs[i]=%.3e\n", ip, xx[ip][1], xx[ip][2], xx[ip][3], funcs[i])
                        #    @printf("normal_ip = [%.3e, %.3e, %.3e]\n", normal_ip[1], normal_ip[2], normal_ip[3])
                        #    @printf("dN = [%.3e, %.3e, %.3e]\n", dN[1], dN[2], dN[3])
                        #    @printf("normals[%d] = [%.3e, %.3e, %.3e]\n", id, normals[id][1], normals[id][2], normals[id][3])
                        #end
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
	    #solids_at_node_i = collect(contact_solids[i]) # to convert from set to array to access using index

	    # skip nodes receive contribution from 1 single solid,
	    # as they are not contact nodes
            length_solids_at_node_i = length(contact_solids[i])
	    if length_solids_at_node_i == 1 continue end 

            # special treatment of normals here, average normal of solid 1 and solid 2
            # just for 2 deformable solids
            s1        = contact_solids[i][1]
            s2        = contact_solids[i][2]

            #@printf("Contacting solids: %d %d\n", s1, s2)

	    nI_body_1 = nodalNormals[s1][i]
	    nI_body_2 = nodalNormals[s2][i]
            #if i==47488
            #    @printf("nodalNormals[%d][%d] = [%.3e, %.3e, %.3e]\n", s2, i, nodalNormals[s2][i][1], nodalNormals[s2][i][2], nodalNormals[s2][i][3])
            #end

	    xI_body_1 = nodalPos[s1][i] / nodalMass[s1][i]
	    xI_body_2 = nodalPos[s2][i] / nodalMass[s2][i]

	    dxI       = xI_body_1 - xI_body_2


	    # if one solid is rigid, then use normal of it
	    # otherwise, use average

	    if typeof(mats[s1]) <: RigidMaterial
	    	nI  = -nI_body_1
                nI /= norm(nI)

		nodalMomentum_1 = nodalMomentum[s2]
		nodalMass_1     = nodalMass[s2]
		velo1           = nodalMomentum_1[i];                       # body 1 velo				
                velocm          = nodalMomentum_S[i]/nodalMass_S[i];        # system velo
                deltaVe1        = velo1 - velocm;
		D1              = dot(deltaVe1,  nI);			    
		
		if ( D1 > 0. ) 
		    dVxNI    = cross(deltaVe1,  nI);			    
		    C1       = norm(dVxNI)		    
		    omega    = dVxNI / C1
		    muPrime1 = min(fric,C1/D1);
                    nodalMomentum_1[i] = velo1 - D1 * nI
                    #if i==47488
                    #    @printf("nodalMomentum_1[%d] = [%f, %f, %f]\n", i, nodalMomentum_1[i][1], nodalMomentum_1[i][2], nodalMomentum_1[i][3])
                    #end
                    #nodalMomentum_1[i] = velo1 - D1*(  nI + muPrime1*cross(nI,omega) );		        
		    
		    #println("approaching\n")
		else			    	
		    #println("separating\n")
		end
	    elseif typeof(mats[s2]) <: RigidMaterial
	        # first solid is deformable, second one is rigid
	    	nI  = -nI_body_2
                #if i==47488
                #    @printf("nI = [%.3e, %.3e, %.3e]\n", nI[1], nI[2], nI[3])
                #end
	    	nI /= norm(nI)
                #if i==47488
                #    @printf("nI = [%f, %f, %f]\n", nI[1], nI[2], nI[3])
                #end

	    	nodalMomentum_1 = nodalMomentum[s1]			
		velo1           = nodalMomentum_1[i];                       # body 1 velo				
                velocm          = nodalMomentum_S[i]/nodalMass_S[i];        # system velo
                deltaVe1        = velo1 - velocm;			    
		D1              = dot(deltaVe1,  nI);
		if ( D1 > 0. ) 
		    dVxNI    = cross(deltaVe1,  nI);			    
		    C1       = norm(dVxNI)		    
		    omega    = dVxNI / C1
		    muPrime1 = min(fric,C1/D1);
		    
                    nodalMomentum_1[i] = velo1 - D1 * nI
                    #if i==47488
                    #    @printf("velo1 = [%f, %f, %f]\n", velo1[1], velo1[2], velo1[3])
                    #    @printf("velocm = [%f, %f, %f]\n", velocm[1], velocm[2], velocm[3])
                    #    @printf("nI = [%f, %f, %f]\n", nI[1], nI[2], nI[3])
                    #    @printf("s2 nodalMomentum_1[%d] = [%f, %f, %f]\n", i, nodalMomentum_1[i][1], nodalMomentum_1[i][2], nodalMomentum_1[i][3])
                    #end
                    #nodalMomentum_1[i] = velo1 - D1*(  nI + muPrime1*cross(nI,omega) );		        
		    
		    #println("approaching\n")
		else			    	
		    #println("separating\n")
		end
		# two deformable solids in contact at node 'i'    
	    else
	    	nI_common  = nI_body_1 - nI_body_2
	    	nI_common /= norm(nI_common)

	    	dxI        = abs(dot(dxI,nI_common))

		#normals    = [nI_common -nI_common]
		normals    = SMatrix{3,2}( nI_common[1], nI_common[2], nI_common[3],
		    	                   -nI_common[1],-nI_common[2],-nI_common[3])
                #println(nI_common)
		for is = 1 : length_solids_at_node_i
	    	    s               = contact_solids[i][is]
	    	    nodalMomentum_1 = nodalMomentum[s]			
		    velo1           = nodalMomentum_1[i];                       # body 1 velo				
	            velocm          = nodalMomentum_S[i]/nodalMass_S[i];        # system velo
	            deltaVe1        = velo1 - velocm;		
	            nI              = @view normals[:,is]	    
		    D1              = dot(deltaVe1,  nI);			    
		    
		    if ( D1 > 0. && dxI < grid.dy/40 ) 
			dVxNI    = cross(deltaVe1,  nI);			    
			C1       = norm(dVxNI)		    
			omega    = dVxNI / C1
			muPrime1 = min(fric,C1/D1);        

	                nodalMomentum_1[i] = velo1 - D1 * nI
	                #nodalMomentum_1[i] = velo1 - D1*(  nI + muPrime1*cross(nI,omega) );		        
			#println("approaching\n")
		    else			    	
			#println("separating\n")
		    end		
	        end
	    end
	end            # end loop over contact nodes
        boundary_nodes = 0


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
	        dxxp       = zeros(3) #xx[ip]
	        #dup       = du[ip]

	        for i = 1:support
	            in = nearPoints[i]; # index of node 'i'
	            Ni = funcs[i]
	            mI = nMass[in]
	            if (mI > alg.tolerance)
		        #invM       = 1.0 / mI
		        vI         = nMomentum[in] #* invM
		        #vvt        += Ni * vI  => too much dissipation
		        vvp       += Ni * (nMomentum[in] - nMomentum0[in])# * invM
                        Dvvp = Ni * (nMomentum[in] - nMomentum0[in])
		        dxxp       += Ni * vI # * dtime
		        #dup       += Ni * vI # * dtime
	            end
	        end
	        vv[ip]      = vvp
	        xx[ip]      += dxxp * dtime
                du[ip]      += dxxp * dtime
	        #du[ip]      = dup * dtime
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
	    vel    = solid.velocity
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

	    cx = cy = cz = 0.
	    min_lambda = 1e22

	    for ip=1:solid.nodeCount 
	  	fint[ip]  = @SVector [0., 0., 0.]
	  	fbody[ip] = @SVector [0., 0., 0.]

	  	cx        = max(cx,abs(du[ip][1]))
	  	cy        = max(cy,abs(du[ip][2]))
	  	cz        = max(cz,abs(du[ip][3]))
	    end

            #@printf("max particle velo: %f %f %f\n", cx, cy, cz)

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
                min_lambda    = min(min_lambda, minimum(real(eigvals(Array(F[ip]))))) # determine the lowest eigenvalue of F[ip]
		J             = det(F[ip])

		if ( J < 0. )
		    @printf("Troubled particle: %f %f \n", xx[ip][1], xx[ip][2])	
		    plotParticles_3D(output,solids,mats,counter)
		    compute_femp(fixes,t)		
		    closeFile(fixes)		
		    @error("J is negative\n")
		end

                Finv          = inv(F[ip])    # no memory alloc
                L             = Fdot*Finv     # no memory alloc when Fdot is a SMatrix
                #@timeit "2"            D             .= (0.5 * dtime) .* (L .+ L')
                D             = (0.5 * dtime) * (L + L')  # no memory alloc
		strain[ip]   += D #0.5 * (vel_grad + vel_grad' + vel_grad * vel_grad')

		#println(strain[ip])
	   	#@timeit "3" update_stress!(stress[ip],mat,strain[ip],F[ip],J,ip)
	   	stress[ip] = update_stress!(stress[ip],mat,strain[ip],D,F[ip],J,ip,dtime)

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

                #################################
                # determine time step based on particle velocity
	    end # end of particle loop


	    cx = cx/dtime + mat.c_dil 
            cy = cy/dtime + mat.c_dil 
            cz = cz/dtime + mat.c_dil 

            (hx, hy, hz) = compute_characteristic_length!(solid)
            #dtime = dt_factor*min(grid.dx/cx,grid.dy/cy,grid.dz/cz)
            @printf("min_lambda=%.4e\t", min_lambda)
            dtime = dt_factor * min(hx/cx, hy/cy, hz/cz) * min_lambda
	end    # end of solid loop

	fix_Dirichlet_solid(solids,data,dtime)

	compute_fext(solids,funcs_surface, normals_surface, weights_surface, gpCoords_surface,data,t)


	##################################################
	# update position of rigid solids
        ##################################################

	if haskey(data, "rigid_body_velo")
	    rigid_solids = data["rigid_body_velo"]::Array{Tuple{Int64,Function},1}

	    for (s,f) in rigid_solids
	 	solid       = solids[s]
		xx          = solid.pos
                vv          = solid.velocity
		#xc          = solid.centroids
		(vex,vey,vez) = f(t)::Tuple{Float64, Float64, Float64} 

		if (vex,vey,vez) == (0.,0.,0.) continue end
		# update all nodes for visualization only
		@inbounds for ip = 1:solid.nodeCount
                    vv[ip]    = @SVector [vex, vey, vez]
		    xx[ip]   += dtime * @SVector [vex,vey,vez]   # no memory alloc, 
		end
                # update the centroids of surface elems
		#@inbounds for ip = 1:solid.surfCount
		#  xc[ip]   += dtime * @SVector [vex,vey,vez]
		#end
	    end
	end

	if haskey(data, "rigid_body_velo_data")
	    rigid_solids = data["rigid_body_velo_data"]::Array{Tuple{Int64,VelocityData},1}

	    for (s,f) in rigid_solids
	 	solid       = solids[s]
		xx          = solid.pos
                vv          = solid.velocity
		#xc          = solid.centroids
		(vex,vey,vez) = f(t)::Tuple{Float64, Float64, Float64} 

		if (vex,vey,vez) == (0.,0.,0.) continue end
		# update all nodes for visualization only
		@inbounds for ip = 1:solid.nodeCount
                    vv[ip]    = @SVector [vex, vey, vez]
		    xx[ip]   += dtime * @SVector [vex,vey,vez]   # no memory alloc,
                    if (ip==1)
                        @printf("\txx[1] = [%4.3e, %4.3e, %4.3e]\t", xx[ip][1], xx[ip][2], xx[ip][3])
                    end
		end
                # update the centroids of surface elems
		#@inbounds for ip = 1:solid.surfCount
		#  xc[ip]   += dtime * @SVector [vex,vey,vez]
		#end
	    end
	end

if haskey(data, "rigid_body_velo_from_file")
    
    rigid_solids = data["rigid_body_velo_from_file"]::Array{Tuple{Int64,Function},1}

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
