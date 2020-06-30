module Util

using Printf
using StaticArrays
using LinearAlgebra

export computeMaxPrincipleStress, report
export get_cells, get_cell

function computeMaxPrincipleStress(stress::MMatrix{2,2,Float64})
  sigmaxx = stress[1,1]
  sigmayy = stress[2,2]
  sigmaxy = stress[1,2]

  delta   = sqrt( (sigmaxx - sigmayy)^2 + 4. * sigmaxy * sigmaxy )

  return 0.5 * ( sigmaxx + sigmayy + delta )
end

function computeMaxPrincipleStress( stress::MMatrix{3,3,Float64} )

  sigmaxx = stress[1,1];
  sigmayy = stress[2,2];
  sigmazz = stress[3,3];
  sigmaxy = stress[1,2];
  sigmayz = stress[2,3];
  sigmaxz = stress[1,3];

  q  =  ( sigmaxx + sigmayy + sigmazz ) * .3333333333333;
  p1 =  sigmaxy*sigmaxy + sigmaxz*sigmaxz + sigmayz*sigmayz;
  p2 = (sigmaxx - q)^2 + (sigmayy - q)^2 + (sigmazz - q)^2 + 2*p1;
  p  =  sqrt(p2/6);

  # B  = (1/p)*(sigma_bar - q*Identity(3))

  B = SMatrix{3,3}( sigmaxx-q, sigmaxy, sigmaxz, sigmaxy, sigmayy-q, sigmayz,  sigmaxz, sigmayz, sigmazz-q)

  if ( abs(p) < 1e-15 ) 
    B = SMatrix{3,3}(0,0,0,0,0,0,0,0,0);
  else
    B /= p; 
  end

  r  = 0.5*det(B);

  phi = 0.
    
  # In exact arithmetic for a symmetric matrix  -1 <= r <= 1
  # but computation error can leave it slightly outside this range.
  if      (r <= -1) 
      phi = 1.0471975511965976;
  elseif (r >= 1)
      phi = 0;
  else
      phi = acos(r)/3;
  end
    
  # the maximum eigenvalue
  
  return  q + 2*p*cos(phi);
end

function report(grid,solids,dtime)
	num = 0
	[num += solids[s].parCount for s=1:length(solids)]
	@printf("Number of material points: %d \n", num)
  @printf("Number of grid points    : %d \n", grid.nodeCount)
  @printf("Grid cell                : %f %f \n", grid.dx, grid.dy)
  @printf("Initial time step        : %.5E \n\n", dtime)
end

function get_ijk(grid,x)
  deltax = grid.dxI;
  deltay = grid.dyI;
  deltaz = grid.dzI;


  i = min(grid.nodeCountX-1, floor(Int64,(x[1]-grid.xmin) * deltax + 1.))
  j = min(grid.nodeCountY-1, floor(Int64,(x[2]-grid.ymin) * deltay + 1.))
  k = min(grid.nodeCountZ-1, floor(Int64,(x[3]-grid.zmin) * deltaz + 1.))

  if i <= 0 i = 1 end
  if j <= 0 j = 1 end
  if k <= 0 k = 1 end
  return (i,j,k)
end

function get_cell(grid,x)
  deltax = grid.dxI;
  deltay = grid.dyI;
  deltaz = grid.dzI;


  i = floor(Int64,(x[1]-grid.xmin) * deltax + 1.)
  j = floor(Int64,(x[2]-grid.ymin) * deltay + 1.)
  k = floor(Int64,(x[3]-grid.zmin) * deltaz + 1.)

  c = i + ( grid.nodeCountX - 1 ) * ( j - 1) + ( grid.nodeCountX - 1 ) * ( grid.nodeCountY - 1 )  * ( k - 1);

  if c >  ( grid.nodeCountX - 1 ) * ( grid.nodeCountY - 1 ) * ( grid.nodeCountZ - 1 ) 
    pritnln(x)
    error("wrong cell")
  end
 
  return i + ( grid.nodeCountX - 1 ) * ( j - 1) + ( grid.nodeCountX - 1 ) * ( grid.nodeCountY - 1 )  * ( k - 1);
end

function get_cells(rad,center,grid)
  # bounding box of the circle
  xmin = center[1] - rad;
  xmax = center[1] + rad;
  ymin = center[2] - rad;
  ymax = center[2] + rad;
  zmin = center[3] - rad;
  zmax = center[3] + rad;

  # indices of cells containing xmin, xmax, ymin, ymax
  (imin,jmin,kmin) = get_ijk(grid,@SVector[xmin,ymin,zmin]);
  (imax,jmax,kmax) = get_ijk(grid,@SVector[xmax,ymax,zmax]);

  cells     = Vector{Int64}(undef,0)
  for k = kmin:kmax
    for j = jmin:jmax
      for i = imin:imax
        id  = i + ( grid.nodeCountX - 1 ) * ( j - 1) + ( grid.nodeCountX - 1 ) * ( grid.nodeCountY - 1 )  * ( k - 1);
        push!(cells,id)
      end
    end
  end

  return cells
end

end

