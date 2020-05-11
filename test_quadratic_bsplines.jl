import PyPlot

pntCount    = 100
node_coords = [0., 1., 2., 3., 4.]      # node coords
points      = LinRange(0.,4.,pntCount)  # points where basis are evaluated
dx          = 1.

function quad_bspline_type1(r)
    if r >= -1.5 && r <= -0.5
        return 0.5*r*r + 1.5*r + 1.125
    elseif r <= 0.
        return 1. + r
    elseif r <= 0.5
        return 1. - r
    elseif r <= 1.5
        return 0.5*r*r -1.5*r + 1.125
    else
        return 0.
    end
end

function quad_bspline_type2(r)
    if r >= -1. && r <= -0.5
        return 1. + r
    elseif r <= 0.5
        return  -r*r +0.75
    elseif r <= 1.5
        return 0.5*r*r - 1.5*r + 1.125
    else
        return 0.
    end
end

function quad_bspline_type3(r)
    if r >= -1.5 && r <= -0.5
        return 0.5*r*r + 1.5*r + 1.125
    elseif r <= 0.5
        return  -r*r + 0.75
    elseif r <= 1.5
        return 0.5*r*r - 1.5*r + 1.125
    else
        return 0.
    end
end


function quad_bspline_type4(r)
    if r >= -1.5 && r <= -0.5
        return 0.5*r*r + 1.5*r + 1.125
    elseif r <= 0.5
        return  -r*r +0.75
    elseif r <= 1.
        return 1. - r
    else
        return 0.
    end
end

basis = zeros(length(points),5)


for p=1:length(points)
    r1         = (points[p] - node_coords[1])/dx
    r2         = (points[p] - node_coords[2])/dx
    r3         = (points[p] - node_coords[3])/dx
    r4         = (points[p] - node_coords[4])/dx
    r5         = (points[p] - node_coords[5])/dx
    basis[p,1] = quad_bspline_type1(r1)
    basis[p,2] = quad_bspline_type2(r2)
    basis[p,3] = quad_bspline_type3(r3)
    basis[p,4] = quad_bspline_type4(r4)
    basis[p,5] = quad_bspline_type1(r5)
end


pyFig_RealTime = PyPlot.figure("MPM 2Disk FinalPlot", figsize=(8/2.54, 4/2.54))
PyPlot.clf()
pyPlot01 = PyPlot.gca()
PyPlot.subplots_adjust(left=0.15, bottom=0.25, right=0.65)
pyPlot01[:grid](b=true, which="both", color="gray", linestyle="-", linewidth=0.5)
pyPlot01[:set_axisbelow](true)
pyPlot01[:set_xlim](0.0, 4.0)
pyPlot01[:set_ylim](0.0, 1.0)
pyPlot01[:set_xlabel]("\$x\$", fontsize=8)
pyPlot01[:set_ylabel]("\$N_I(x)\$", fontsize=8)
pyPlot01[:set_xticks](collect(0.0:1.0:4.0))
pyPlot01[:tick_params](axis="both", which="major", labelsize=8)
pyPlot01[:set_yticks](collect(0.0:0.2:1.0))
PyPlot.plot(points, basis[:,1], "-",  c="red", label="\$ K \$", linewidth=1.0)
PyPlot.plot(points, basis[:,2], "-",  c="black", label="\$ K \$", linewidth=1.0)
PyPlot.plot(points, basis[:,3], "-",  c="blue", label="\$ K \$", linewidth=1.0)
PyPlot.plot(points, basis[:,4], "-",  c="black", label="\$ K \$", linewidth=1.0)
PyPlot.plot(points, basis[:,5], "-",  c="red", label="\$ K \$", linewidth=1.0)
#PyPlot.hold(true)
PyPlot.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=8)
