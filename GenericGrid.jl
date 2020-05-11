struct Grid2D
    lx        :: Float64              # length in x dir
    ly        :: Float64              # length in y dir
    dx        :: Float64              # cell size in y dir
    dy        :: Float64              # cell size in y dir
    dxI       :: Float64              # inverse of cell size in x dir
    dyI       :: Float64              # inverse of cell size in y dir
    nodeCount :: Int64                # number of grid nodes
    nodeCountX:: Int64                # number of grid nodes along x dir
    nodeCountY:: Int64                # number of grid nodes along y dir

    mass      :: Vector{Float64}
    pos       :: Vector{SVector{2,Float64}}
    momentum  :: Vector{MVector{2,Float64}}
    momentum2 :: Vector{MVector{2,Float64}}  # for MUSL algorithm
    force     :: Vector{MVector{2,Float64}}

    fixedXNodes::Vector{Int64}         # 1D indices of all nodes fixed in X dir.: 1 = fixed, 0: free
    fixedYNodes::Vector{Int64}         # 1D indices of all nodes fixed in Y dir.

    # constructor, GL_x is length of the grid in x dir
    # iN_x: number of nodes in x dir
    function Grid2D(fGL_x, fGL_y, iN_x, iN_y)
       dx  = fGL_x / Float64(iN_x - 1.0)
       dy  = fGL_y / Float64(iN_y - 1.0)
       dxI = 1.0 / dx
       dyI = 1.0 / dy


       mass         = fill(0,iN_x*iN_y)
       momentum     = fill(zeros(2),iN_x*iN_y)
       momentum2    = fill(zeros(2),iN_x*iN_y)
       force        = fill(zeros(2),iN_x*iN_y)
       pos          = fill(zeros(2),iN_x*iN_y)

       idx          = fill(0,iN_x*iN_y)
       idy          = fill(0,iN_x*iN_y)

       for j=1:iN_y
          for i=1:iN_x
             x     = (i-1) * dx
             y     = (j-1) * dy
             index = index2DTo1D(i, j, iN_x, iN_y)
             pos[index] = [x, y];
          end
      end
      new(fGL_x, fGL_y, dx, dy, dxI, dyI, iN_x*iN_y, iN_x, iN_y, mass, pos, momentum, momentum2, force,idx,idy)
    end
 end
# ----------------------------------------------------------------------
