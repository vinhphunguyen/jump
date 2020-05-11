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
