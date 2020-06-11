module BodyForce

using StaticArrays

abstract type BodyForceType end

struct XDirection end
struct YDirection end
struct ZDirection end

struct ConstantBodyForce1D <: BodyForceType
    g::Float64
end

# Body force for the 1D MMS
struct AxisAlignBodyForce1D  <: BodyForceType
    G::Float64
    E::Float64
    nu::Float64
    rho::Float64
    c::Float64
    lambda::Float64
    mu    ::Float64
    function AxisAlignBodyForce1D(G,E,nu,rho)
        lambda = E*nu/(1+nu)/(1-2*nu)
        mu     = E/2/(1+nu)
        c      = sqrt(E/rho)
        new(G,E,nu,rho,c,lambda,mu)
    end
end

struct VortexBodyForce2D  <: BodyForceType
    G::Float64
    T::Float64
    Ri::Float64
    Ro::Float64
    rho0::Float64
    mu    ::Float64
    R    ::Float64
    Rio2 ::Float64
    Rio4 ::Float64
    RimRo ::Float64
    xcenter ::Float64
    ycenter ::Float64
    function VortexBodyForce2D(G,T,Ri,Ro,rho0,mu,xcenter,ycenter)
        R      = (Ri+Ro)/2
        Rio2   = (Ri-Ro)^2
        Rio4   = (Ri-Ro)^4
        RimRo  = (Ri-Ro)
        new(G,T,Ri,Ro,rho0,mu,R,Rio2,Rio4,RimRo,xcenter,ycenter)
    end
end



struct ConstantBodyForce2D <: BodyForceType
    g::SVector{2}
end

struct ConstantBodyForce3D <: BodyForceType
    g::SVector{3}
end

function (b::ConstantBodyForce1D)(x,t)
    return -b.g*rho
end

function (b::AxisAlignBodyForce1D)(x,t)
    cpit = b.c*pi*t
    u    = b.G*sin(pi*x)*sin(cpit)
    F    = 1.0 + pi*b.G*cos(pi*x)*sin(cpit)
    F2i  = 1/F/F
    body = pi*pi*u/b.rho*(b.lambda*F2i*(1.0 - log(F)) + b.mu*(1.0 + F2i) - b.E )
    return body#*b.rho
end

function evaluate_stress(b::AxisAlignBodyForce1D,x,t)
    cpit = b.c*pi*t
    F    = 1.0 + pi*b.G*cos(pi*x)*sin(cpit)
    F2i  = 1/F
    stress = F2i * ( b.lambda*log(F) + b.mu * (F*F-1.0) )
    return stress
end

function (b::ConstantBodyForce2D)(body,x,t)
     body .= b.g
end

function (b::ConstantBodyForce3D)(body,x,t)
    body  .= b.g
end


function (b::VortexBodyForce2D)(body,xp,t)
    G     = b.G
    T     = b.T
    R     = b.R
    Ri    = b.Ri
    Ro    = b.Ro

    x     = xp[1]-b.xcenter
    y     = xp[2]-b.ycenter
    r     = sqrt(x*x+y*y)
    theta = atan(y,x)

    R     = (Ri+Ro)/2
    s     = (r - R)/(Ri-Ro)
    h     = 1 - 8*((r - R)/(Ri-Ro))^2 +16*((r - R)/(Ri-Ro))^4
    hp    = - 16*(r - R)/(Ri-Ro)^2 + 16*4*(r - R)^3/(Ri-Ro)^4
    hpp   = - 16/(Ri-Ro)^2 + 16*4*3*(r - R)^2/(Ri-Ro)^4
    g     = G * sin(π*t/T)
    gp    = G*π/T*cos(π *t/T)
    gpp   = -π^2/T^2*g
    alpha = g*h
    mdr   = b.mu/b.rho0

    br    = ( mdr*(3*g*hp+r*g*hpp) - r*gpp*h)*sin(alpha) + (mdr*r*(g*hp)^2 - r*(gp*h)^2)*cos(alpha)
    bt    = (-mdr*(3*g*hp+r*g*hpp) + r*gpp*h)*cos(alpha) + (mdr*r*(g*hp)^2 + r*(gp*h)^2)*sin(alpha)

    body .= @SVector[br*cos(theta) - bt*sin(theta), br*sin(theta) + bt*cos(theta)]
end    

export BodyForceType, ConstantBodyForce1D, ConstantBodyForce3D, AxisAlignBodyForce1D, VortexBodyForce2D,
       ConstantBodyForce2D, evaluate_stress

end
