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

module Grid
  using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")
  using Printf

  # ----------------------------------------------------------------------
  struct Grid1D
     lx        :: Float64              # length
     dx        :: Float64              # cell size in y dir
     dxI       :: Float64              # inverse of cell size in x dir
     nodeCount :: Int64                # number of grid nodes

     mass      :: Vector{Float64}
     pos       :: Vector{Float64}
     pos0      :: Vector{Float64}
     momentum  :: Vector{Float64}
     momentum2 :: Vector{Float64}  # for MUSL algorithm
     force     :: Vector{Float64}

     fixedNodes::Vector{Int64}         # 1D indices of all nodes fixed in X dir.: 1 = fixed, 0: free

     # constructor, length and number of nodes
     function Grid1D(fGL_x, iN_x)
        dx           = fGL_x / Float64(iN_x - 1.0)
        dxI          = 1.0 / dx

        mass         = fill(0.,iN_x)
        momentum     = fill(0.,iN_x)
        momentum2    = fill(0.,iN_x)
        force        = fill(0.,iN_x)
        pos          = fill(0.,iN_x)
        pos0         = fill(0.,iN_x)
        idx          = fill(0,iN_x)

        for i=1:iN_x
           x      = (i-1) * dx
           pos[i] =  x;
           pos0[i]=  x;
        end

       return (fGL_x, dx,dxI, iN_x, mass, pos, pos0, momentum,
                     momentum2, force,idx)
     end
  end
  # ----------------------------------------------------------------------

  # ----------------------------------------------------------------------
  #grid container
   struct Grid2D
      xmin      ::Float64               # for TLMPM, grid origin is not always (0,0)
      ymin      ::Float64
      lx        :: Float64              # length in x dir
      ly        :: Float64              # length in y dir
      dx        :: Float64              # cell size in y dir
      dy        :: Float64              # cell size in y dir
      dxI       :: Float64              # inverse of cell size in x dir
      dyI       :: Float64              # inverse of cell size in y dir
      nodeCount :: Int64                # number of grid nodes
      nodeCountX:: Int64                # number of grid nodes along x dir
      nodeCountY:: Int64                # number of grid nodes along y dir
      cellCount :: Int64                # number of grid cells

      mass      :: Vector{Float64}
      pos       :: Vector{SVector{2,Float64}}
      momentum0 :: Vector{SVector{2,Float64}}  # momentum at the beginning of the step
      momentum  :: Vector{SVector{2,Float64}}  # momentum at the end of the step (tilde)
      momentum2 :: Vector{SVector{2,Float64}}  # for MUSL algorithm
      force     :: Vector{SVector{2,Float64}}

      fixedXNodes::Vector{Int64}         # 1D indices of all nodes fixed in X dir.: 1 = fixed, 0: free
      fixedYNodes::Vector{Int64}         # 1D indices of all nodes fixed in Y dir.

      damage     :: Vector{Float64}      # nodal damage for phase field
      damage0    :: Vector{Float64}      # nodal damage for phase field
      pfForce    :: Vector{Float64}      # nodal force for phase field

      # constructor, GL_x is length of the grid in x dir
      # iN_x: number of nodes in x dir
      function Grid2D(xmin, xmax, ymin, ymax, iN_x, iN_y)
         dx  = (xmax-xmin) / Float64(iN_x - 1.0)
         dy  = (ymax-ymin) / Float64(iN_y - 1.0)
         dxI = 1.0 / dx
         dyI = 1.0 / dy

         nodeCount    = iN_x*iN_y
         mass         = fill(0.,nodeCount)
         momentum0    = [@SVector [0., 0.] for _ in 1:nodeCount]
         momentum     = [@SVector [0., 0.] for _ in 1:nodeCount]
         momentum2    = [@SVector [0., 0.] for _ in 1:nodeCount]
         force        = [@SVector [0., 0.] for _ in 1:nodeCount]
         pos          = [@SVector [0., 0.] for _ in 1:nodeCount]

         idx          = fill(0,nodeCount)
         idy          = fill(0,nodeCount)

         for j=1:iN_y
            for i=1:iN_x
               x     = (i-1) * dx
               y     = (j-1) * dy
               index = index2DTo1D(i, j, iN_x, iN_y)
               pos[index] = @SVector [xmin+x, ymin+y]
            end
        end
        return new(xmin, ymin, xmax-xmin, ymax-ymin,dx, dy, dxI, dyI, nodeCount,
            iN_x, iN_y, (iN_x-1)*(iN_y-1) ,mass, pos, momentum0,
            momentum, momentum2, force,idx,idy,copy(mass),copy(mass),copy(mass))
      end
   end
  # ----------------------------------------------------------------------

  struct Grid3D
     xmin      :: Float64
     ymin      :: Float64
     zmin      :: Float64
     lx        :: Float64              # length in x dir
     ly        :: Float64              # length in y dir
     lz        :: Float64              # length in z dir
     dx        :: Float64              # cell size in y dir
     dy        :: Float64              # cell size in y dir
     dz        :: Float64              # cell size in z dir
     dxI       :: Float64              # inverse of cell size in x dir
     dyI       :: Float64              # inverse of cell size in y dir
     dzI       :: Float64              # inverse of cell size in z dir
     nodeCount :: Int64                # number of grid nodes
     nodeCountX:: Int64                # number of grid nodes along x dir
     nodeCountY:: Int64                # number of grid nodes along y dir
     nodeCountZ:: Int64                # number of grid nodes along  dir
     nodeCountXY:: Int64               # number of grid nodes one XY plane

     mass      :: Vector{Float64}
     pos       :: Vector{SVector{3,Float64}}
     momentum0 :: Vector{MVector{3,Float64}}
     momentum  :: Vector{MVector{3,Float64}}
     momentum2 :: Vector{MVector{3,Float64}}  # for MUSL algorithm
     force     :: Vector{MVector{3,Float64}}

     fixedNodes :: Array{Int64,2}
     bottomNodes:: Vector{Int64}
     topNodes   :: Vector{Int64}
     leftNodes  :: Vector{Int64}
     rightNodes :: Vector{Int64}
     frontNodes :: Vector{Int64}
     backNodes  :: Vector{Int64}

     # constructor, GL_x is length of the grid in x dir
     # iN_x: number of nodes in x dir
     function Grid3D(xmin,xmax,ymin,ymax,zmin,zmax, iN_x, iN_y, iN_z)
        dx  = (xmax-xmin) / Float64(iN_x - 1.0)
        dy  = (ymax-ymin) / Float64(iN_y - 1.0)
        dz  = (zmax-zmin) / Float64(iN_z - 1.0)
        dxI = 1.0 / dx
        dyI = 1.0 / dy
        dzI = 1.0 / dz

        nodeCount    = iN_x*iN_y*iN_z

        mass         = fill(0,nodeCount)
        momentum0    = fill(zeros(3),nodeCount)
        momentum     = fill(zeros(3),nodeCount)
        momentum2    = fill(zeros(3),nodeCount)
        force        = fill(zeros(3),nodeCount)
        pos          = fill(zeros(3),nodeCount)

        bNodes       = Vector{Int64}(undef,iN_x*iN_z)
        tNodes       = Vector{Int64}(undef,iN_x*iN_z)
        lNodes       = Vector{Int64}(undef,iN_y*iN_z)
        rNodes       = Vector{Int64}(undef,iN_y*iN_z)
        fNodes       = Vector{Int64}(undef,iN_x*iN_y)
        baNodes      = Vector{Int64}(undef,iN_x*iN_y)
 

        fixedNodes   = Array{Int64}(undef, 3, nodeCount)
        fixedNodes  .= 0

        c = 1
        for k=1:iN_z
          @inbounds for i=1:iN_x
            ii         = i                 + (k-1)*iN_x*iN_y 
            jj         = i + iN_x*(iN_y-1) + (k-1)*iN_x*iN_y              
            bNodes[c]  = ii  
            tNodes[c]  = jj
            c         += 1    
          end
        end

        c = 1
        for k=1:iN_z
          for i=1:iN_y
            rNodes[c] = iN_x*i       + (k-1)*iN_x*iN_y  
            lNodes[c] = iN_x*(i-1)+1 + (k-1)*iN_x*iN_y 
            c         += 1    
          end
        end

       @inbounds for i=1:iN_x*iN_y  
          fNodes[i]  = i
          baNodes[i] = i+(iN_z-1)*iN_x*iN_y 
        end
  

        # build grid node coords
        for k=1:iN_z
           for j=1:iN_y
              for i=1:iN_x
                 x     = (i-1) * dx
                 y     = (j-1) * dy
                 z     = (k-1) * dz
                 index = index3DTo1D(i, j, k, iN_x, iN_y, iN_z)

                 pos[index] = @SVector [xmin+x, ymin+y, zmin+z]
              end
           end
      end
      new(xmin, ymin, zmin, xmax-xmin, ymax-ymin, zmax-zmin, dx, dy, dz, dxI, dyI, dzI, nodeCount, iN_x,
          iN_y,iN_z, iN_x*iN_y, mass, pos, momentum0,momentum, momentum2, force,fixedNodes,bNodes,tNodes,lNodes,rNodes,fNodes,baNodes)
     end
  end
# ----------------------------------------------------------------------

  function index2DTo1D(i::Int64, j::Int64, nColumns::Int64, nRows::Int64)
    index = nColumns*(j-1) + i

    if(index > nRows*nColumns || index < 1)
      @printf("Index out of bounds: i, j: %d, %d \n", i, j)
    end

    return(Int64(index))
  end

  function index3DTo1D(i, j, k, nColumns, nRows, nLayers)
     index = nColumns*nRows*(k-1) + nColumns*(j-1) + i

     if(index > nRows*nColumns*nLayers)
        @printf("Index out of bounds")
     end

     return(Int(index))
  end

  function fixXForBottom(grid::Grid2D)
    xnodes = grid.fixedXNodes
    @inbounds for i=1:grid.nodeCountX
      xnodes[i] = 1
    end
  end

  function fixYForBottom(grid::Grid2D)
    ynodes = grid.fixedYNodes
    @inbounds for i=1:grid.nodeCountX
      ynodes[i] = 1
    end
  end


  function fixXForTop(grid::Grid2D;ghostcell=false)
    xnodes = grid.fixedXNodes
    @inbounds for i=1:grid.nodeCountX
      xnodes[i+grid.nodeCountX*(grid.nodeCountY-1)] = 1
    end
    if (ghostcell)
         @inbounds for i=1:grid.nodeCountX
            xnodes[i+grid.nodeCountX*(grid.nodeCountY-1)-grid.nodeCountX] = 1
         end
    end
  end

  function fixYForTop(grid::Grid2D;ghostcell=false)
    ynodes = grid.fixedYNodes
    @inbounds for i=1:grid.nodeCountX
      ynodes[i+grid.nodeCountX*(grid.nodeCountY-1)] = 1
    end
    if (ghostcell)
        @inbounds for i=1:grid.nodeCountX
           ynodes[i+grid.nodeCountX*(grid.nodeCountY-1)-grid.nodeCountX] = 1
        end
    end
  end

  function fixXForLeft(grid::Grid2D;ghostcell=false)
    xnodes = grid.fixedXNodes
    @inbounds for i=1:grid.nodeCountY
      xnodes[grid.nodeCountX*(i-1)+1] = 1
    end
    if (ghostcell)
         @inbounds for i=1:grid.nodeCountY
            xnodes[grid.nodeCountX*(i-1)+2] = 1
         end
    end
  end

  function fixYForLeft(grid::Grid2D;ghostcell=false)
    ynodes = grid.fixedYNodes
    @inbounds for i=1:grid.nodeCountY
      ynodes[grid.nodeCountX*(i-1)+1] = 1
    end
    if (ghostcell)
         @inbounds for i=1:grid.nodeCountY
            ynodes[grid.nodeCountX*(i-1)+2] = 1
         end
    end
  end

  function fixXForRight(grid::Grid2D)
    xnodes = grid.fixedXNodes
    for i=1:grid.nodeCountY
      xnodes[grid.nodeCountX*i] = 1
    end
  end

  function fixYForRight(grid::Grid2D)
    ynodes = grid.fixedYNodes
    @inbounds for i=1:grid.nodeCountY
      ynodes[grid.nodeCountX*i] = 1
    end
  end

  function fixXForLeft(grid::Grid1D;ghostcell=false)
    xnodes    = grid.fixedNodes
    xnodes[1] = 1
  end

  function fixXForRight(grid::Grid1D;ghostcell=false)
    xnodes                 = grid.fixedNodes
    xnodes[grid.nodeCount] = 1
  end

  

  # find the neighbors of a given element "elemId"
  # in a structured grid of numx x numy elements
  function getNeighbors(elemId, grid::Grid2D)
      numx = grid.nodeCountX-1
      numy = grid.nodeCountY-1
      row  = ceil(elemId/numx)
      col  = rem(elemId,numx)

      if col == 0
          col = numx
      end

      if     row  == 1
          rows = [1 2]
      elseif row == numy
          rows = [numy-1 numy]
      else
          rows = [row-1 row row + 1]
      end

      if     col == 1
          cols = [1 2];
      elseif col == numx
          cols = [numx-1 numx]
      else
          cols = [col-1 col col + 1]
      end
      neighbors= repeat(0:0, inner=length(rows)*length(cols))
      ii = 1
      for i=1:length(rows)
        id = rows[i]
        for j=1:length(cols)
          jd = cols[j]
          neighbors[ii] = (id-1)*numx+jd
          ii+=1
        end
      end
      return neighbors
  end

  export Grid1D, Grid2D, Grid3D
  export index2DTo1D, index3DTo1D, getNeighbors,
         fixXForBottom, fixYForBottom, fixZForBottom, fixForBottom, fixXForTop, fixYForTop, fixXForLeft, fixYForLeft, fixZForLeft, fixXForRight, fixYForRight, fixZForRight, fixForTop

end
