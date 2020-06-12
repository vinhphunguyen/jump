module Util

using Printf

export computeMaxPrincipleStress2D, report

function computeMaxPrincipleStress2D(stress)
  sigmaxx = stress[1,1]
  sigmayy = stress[2,2]
  sigmaxy = stress[1,2]

  delta   = sqrt( (sigmaxx - sigmayy)^2 + 4. * sigmaxy * sigmaxy )

  return 0.5 * ( sigmaxx + sigmayy + delta )
end

function report(grid,solids,dtime)
	num = 0
	[num += solids[s].parCount for s=1:length(solids)]
	@printf("Number of material points: %d \n", num)
  @printf("Number of grid points    : %d \n", grid.nodeCount)
  @printf("Grid cell                : %f %f \n", grid.dx, grid.dy)
  @printf("Initial time step        : %.5E \n\n", dtime)
end


end
