import PyPlot
using  Solid

function plotParticles(solids::Vector{Solid2D},fileName)
	array_x         = [toXArray(solids[i]) for i = 1:length(solids)]
	array_y         = [toYArray(solids[i]) for i = 1:length(solids)]
	#    array_color     = Array{Real}(iMaterialPoints, 3)
	# array_size      = Array{Real}(iMaterialPoints, 1)
	#    for iIndex in 1:1:iMaterialPoints
	#      array_color[iIndex, :] = [1.0, 0.0, 0.0]#[thisGrid.GridPoints[iIndex].fMass/iMaterialPoints, 0.0, 0.0]
	#   array_size[iIndex, :] = [5.0]
	#    end

	pyPlot01 = PyPlot.gca()
	# pyPlot01 = PyPlot.subplot2grid((1,1), (0,0), colspan=1, rowspan=1, aspect="equal")
	PyPlot.scatter(array_x, array_y, lw=0)
	pyPlot01[:spines]["top"][:set_color]("gray")
	pyPlot01[:spines]["right"][:set_color]("gray")
	pyPlot01[:spines]["bottom"][:set_color]("gray")
	pyPlot01[:spines]["left"][:set_color]("gray")
	# pyPlot01[:axhline](linewidth=4, color="g")
	# pyPlot01[:axvline](linewidth=4, color="g")
	pyPlot01[:set_xlim](0.0, 1.0)
	pyPlot01[:set_ylim](0.0, 1.0)
	# pyPlot01[:set_xlabel]("")
	# pyPlot01[:set_ylabel]("")
	pyPlot01[:grid](b=true, which="both", color="gray", linestyle="-", linewidth=0.5)
	pyPlot01[:set_axisbelow](true)
	pyPlot01[:set_xticks]([])# empty to have no major ticks and grids
	pyPlot01[:set_xticks](collect(0.0:0.05:1.0),minor=true)
	pyPlot01[:set_yticks]([])# empty to have no major ticks and grids
	pyPlot01[:set_yticks](collect(0.0:0.05:1.0),minor=true)

	PyPlot.show()
	# PyPlot.hold(true)

	strFileName = "image.png"
	PyPlot.savefig(fileName, bbox_inches="tight")
end