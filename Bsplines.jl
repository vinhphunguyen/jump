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


#-----------------------------------------------------------------
# Modified cubic Bsplines
#-----------------------------------------------------------------

function quad_bspline_type1(r)
    if r < -1.5 || r > 1.5 return (0.,0.) end

    if r >= -1.5 && r <= -0.5
        N  =  0.5*r*r + 1.5*r + 1.125
        dN =  r + 1.5
    elseif r <= 0.
        N = 1. + r
        dN = 1.
    elseif r <= 0.5
        N  = 1. - r
        dN = -1.
    else#if r <= 1.5
        N  = 0.5*r*r -1.5*r + 1.125
        dN = r -1.5
    end
    return (N,dN)
end

function quad_bspline_type2(r)
    if r < -1.0 || r > 1.5 return (0.,0.) end

    if r >= -1. && r <= -0.5
        N  = 1. + r
        dN = 1.
    elseif r <= 0.5
        N  =  -r*r +0.75
        dN =  -2*r
    else#if r <= 1.5
        N  = 0.5*r*r - 1.5*r + 1.125
        dN = r - 1.5
    end
    return (N,dN)
end

function quad_bspline_type3(r)
    if r < -1.5 || r > 1.5 return (0.,0.) end

    if r >= -1.5 && r <= -0.5
        N  = 0.5*r*r + 1.5*r + 1.125
        dN = r + 1.5
    elseif r <= 0.5
        N  =  -r*r + 0.75
        dN =  -2*r
    else#if r <= 1.5
        N  = 0.5*r*r - 1.5*r + 1.125
        dN = r - 1.5
    end
    return (N,dN)
end


function quad_bspline_type4(r)
    if r < -1.5 || r > 1. return (0.,0.) end
    if r >= -1.5 && r <= -0.5
        N  = 0.5*r*r + 1.5*r + 1.125
        dN = r + 1.5
    elseif r <= 0.5
        N  =  -r*r +0.75
        dN =  -2*r
    else#if r <= 1.
        N  = 1. - r
        dN = -1.
    end
    return (N,dN)
end

#-----------------------------------------------------------------
# Modified cubic Bsplines
#-----------------------------------------------------------------

function cubic_bspline_type1(r)
    if r < -2. || r > 2. return (0.,0.) end

    if r >= -2.0 && r <= -1.0
        r2 = r*r 
        N  =  0.16666666666666666*r*r2 + r2 + 2*r + 1.3333333333333333
        dN =  0.5*r2 + 2*r + 2.
    elseif -1 <= r && r <= 0.
        r2 = r*r
        N  = -0.16666666666666666*r*r2 + r + 1
        dN = -0.5*r2 + 1.
    elseif 0 <= r && r <= 1.0
        r2 = r*r
        N  = 0.16666666666666666*r*r2 - r + 1
        dN = 0.5*r2 - 1.0
    else#if r <= 1.5
        r2 = r*r
        N  = -0.16666666666666666*r*r2 + r2 -2*r + 1.3333333333333333
        dN = -0.5*r2 + 2*r - 2.
    end
    return (N,dN)
end


function cubic_bspline_type2(r)
    if r < -1. || r > 2. return (0.,0.) end

    if r >= -1.0 && r <= 0.
        r2 = r*r 
        N  =  -0.3333333333333333*r*r2 - r2 +  0.6666666666666666
        dN =  -r2 - 2*r 
    elseif 0. <= r && r <= 1.
        r2 = r*r
        N  = 0.5*r*r2 - r2 + 0.6666666666666666
        dN = 1.5*r2 -2*r
    elseif 1 <= r && r <= 2.0
        r2 = r*r
         N  = -0.16666666666666666*r*r2 + r2 -2*r + 1.3333333333333333
        dN = -0.5*r2 + 2*r - 2
    end
    return (N,dN)
end

function cubic_bspline_type3(r)
    if r < -2. || r > 2. return (0.,0.) end

    if r >= -2.0 && r <= -1.0
        r2 = r*r 
        N  =  0.16666666666666666*r*r2 + r2 + 2*r + 1.3333333333333333
        dN =  0.5*r2 + 2*r + 2
    elseif -1. <= r && r <= 0.
        r2 = r*r
        N  = -0.5*r*r2 - r2 + 0.6666666666666666
        dN = -1.5*r2 -2*r
    elseif 0 <= r && r <= 1.0
         r2 = r*r
         N  = 0.5*r*r2 - r2 + 0.6666666666666666
        dN  = 1.5*r2 -2*r
    elseif 1. <= r && r <= 2.0
        r2 = r*r
        N  = -0.16666666666666666*r*r2 + r2 -2*r + 1.3333333333333333
        dN = -0.5*r2 + 2*r - 2      
    end
    return (N,dN)
end

function cubic_bspline_type4(r)
    if r < -2. || r > 1. return (0.,0.) end

    if r >= -2.0 && r <= -1.0
        r2 = r*r 
        N  = 0.16666666666666666*r*r2 + r2 +2*r + 1.3333333333333333
        dN = 0.5*r2 + 2*r + 2   
    elseif -1. <= r && r <= 0.
        r2 = r*r
        N  = -0.5*r*r2 - r2 + 0.6666666666666666
        dN = -1.5*r2 -2*r
    elseif 0 <= r && r <= 1.0
        r2 = r*r
        N  =  0.3333333333333333*r*r2 - r2 +  0.6666666667
        dN =  r2 - 2*r 
    end
    return (N,dN)
end
