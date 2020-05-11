module Util

export computeMaxPrincipleStress2D

function computeMaxPrincipleStress2D(stress)
  sigmaxx = stress[1,1]
  sigmayy = stress[2,2]
  sigmaxy = stress[1,2]

  delta   = sqrt( (sigmaxx - sigmayy)^2 + 4. * sigmaxy * sigmaxy )

  return 0.5 * ( sigmaxx + sigmayy + delta )
end

end
