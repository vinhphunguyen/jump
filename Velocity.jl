module Velocity
export VelocityData
######################################################################
# Functor that stores the velocity data and the function
# that interpolates between the data points
######################################################################

struct VelocityData
    T::Array{Float64,1}
    vx::Array{Float64,1}
    vy::Array{Float64,1}
    vz::Array{Float64,1}

    function VelocityData(T, vx, vy, vz)
        if T[1]!=0
            error("VelocityData.T[1] should be equal to 0.")
        end
        if (length(T)!=length(vx)) || (length(T)!=length(vy)) || (length(T)!=length(vz))
            error("VelocityData: all arrays should be of the same length!")
        end
        new(T, vx, vy, vz)
    end
end

# Function that interpolates the velocity between the data points stored in the VelocityData arrays
function (v::VelocityData)(t)
    Idx = findall((v.T.-t).>=0) # Return the indeces of T such that T(i)-t > 0.
    # If t is larger than the end value of T,
    if length(Idx) == 0
        vx = v.vx[end]
        vy = v.vy[end]
        vz = v.vz[end]
    else
        i = Idx[1]
        vx = v.vx[i]
        vy = v.vy[i]
        vz = v.vz[i]
    end
    return (vx,vy,vz)
end

end
