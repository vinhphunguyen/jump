module AbaqusReader

using Nullables

include("parse_mesh.jl")
include("keyword_register.jl")
include("parse_model.jl")
include("create_surface_elements.jl")
#include("abaqus_download.jl")

export abaqus_read_mesh, abaqus_read_model, create_surface_elements
export abaqus_download

end