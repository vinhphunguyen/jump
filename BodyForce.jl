module BodyForce

abstract type BodyForceType end

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

struct ConstantBodyForce2D <: BodyForceType
    g::Float64
end

struct ConstantBodyForce3D <: BodyForceType
    g::Float64
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
    body[2] = -b.g
end

function (b::ConstantBodyForce3D)(body,x,t)
    body[3] = -b.g
end

export BodyForceType, ConstantBodyForce1D, ConstantBodyForce3D, AxisAlignBodyForce1D,
       ConstantBodyForce2D, evaluate_stress

end
