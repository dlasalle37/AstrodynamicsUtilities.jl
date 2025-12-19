abstract type AbstractStateRepresentation end


# Containers
struct CartesianState{T<:Union{Nothing,AbstractVector}} <: AbstractStateRepresentation
    vec::T

    function CartesianState(vec::TP) where TP
        return new{TP}(vec)
    end

    function CartesianState()
        return new{Nothing}(nothing)
    end

end

struct KeplerianState{T<:Union{Nothing,AbstractVector}} <: AbstractStateRepresentation
    vec::T
    
    function KeplerianState(vec::TP) where TP
        return new{TP}(vec)
    end

    function KeplerianState()
        return new{Nothing}(nothing)
    end
end

# Exported Functions

"""
    function transform(old::CartesianState{TP}, new::KeplerianState{Nothing}) where TP<:AbstractVector
Transforms an inertial Cartesian state to an inertial Keplerian state. Uses the 'rv2coe' algorithm from Fundamentals of Astro. Vallado (2001).

# Inputs:
    old: Cartesian state vector container
    new_type: Flag to define what type to transform to. This method is to transform to KeplerianState().
    mu: gravitational parameter of central body (e.g., 398600.4418) for Earth in km^3/s^2

# Returns:
    KeplerianState(vec): Keplerian container with transformed vec
"""
function transform(old::CartesianState{TP}, new_type::KeplerianState{Nothing}, mu) where TP<:AbstractVector

    
    vec = rv2coe(old.vec, mu)


    return KeplerianState(vec)
end

"""
    function transform(old::KeplerianState{TP}, new_type::CartesianState{Nothing}, mu) where TP<:AbstractVector
Transforms an inertial keplerian state to an inertial ceplerian state. Uses the 'coe2rv' algorithm from Fundamentals of Astro. Vallado (2001).

# Inputs:
    old: keplerian state vector container
    new_type: Flag to define what type to transform to. This method is to transform to CartesianState().
    mu: gravitational parameter of central body (e.g., 398600.4418) for Earth in km^3/s^2

# Returns:
    KeplerianState(vec): Keplerian container with transformed vec
"""
function transform(old::KeplerianState{TP}, new_type::CartesianState{Nothing}, mu) where TP<:AbstractVector

    
    vec = coe2rv(old.vec, mu)


    return CartesianState(vec)
end



# Inner functions
"""
coe2rv: Calculate the inertial position and velocity vectors from Keplerian classical orbital elements.

# Inputs:\n
    `vec`: 6-element vector containing:
        a: semi major axis (km)
        e: eccentricity
        inc: inclination (radians)
        RAAN: right-ascension (radians)
        w: argument of perigee (radians)
        trueAnom: true anomaly (radians)
    `mu`: gravitational parameter of central body (km^3/s^2)
# Returns:\n
    `rv`: 6-element vector of position and velocity
"""
function coe2rv(vec, mu)
    a, e, i, RAAN, w, nu = view(vec, 1:6)
    p = a*(1-e^2); #semi-latus rectum
    r = p / (1+e*cos(nu)); #scalar r from trajectory

    r_p = SA[r*cos(nu); r*sin(nu); 0]; #perifocal position vector


    v_p = SA[-sqrt(mu/p) * sin(nu); sqrt(mu/p)*(e+cos(nu)); 0]; #perifocal v

    #need to transform r_p and v_p  to inertial
    #use p->n transfer matrix, A. (eq126):

    A11 = cos(RAAN)*cos(w) - sin(RAAN)*sin(w)*cos(i);
    A12 = -cos(RAAN)*sin(w) - sin(RAAN)* cos(w)*cos(i);
    A13 = sin(RAAN)*sin(i);
    A21 = sin(RAAN)*cos(w) + cos(RAAN)*sin(w)*cos(i);
    A22 = -sin(RAAN)*sin(w) + cos(RAAN)*cos(w)*cos(i);
    A23 = -cos(RAAN)*sin(i);
    A31 = sin(w)*sin(i);
    A32 = cos(w)*sin(i);
    A33 = cos(i);

    A = SA[A11 A12 A13;
        A21 A22 A23;
        A31 A32 A33];
    #inertial vectors:

    r = A*r_p; 
    v = A*v_p;
    return SA[r...; v...]
end

"""
    rv2coe(r, v, mu)

Calculate the 6 classical keplerian orbital elements from r, v, mu (under the two body problem)

# Inputs:\n
    r: inertial position vector (km)
    v: inertial velocity vector (km/s)
    mu: gravitational parameter of central body (km^3/s^2)

# returns:\n
    coe: vector of 6 Keplerian elements, which are:\n
        a: semi major axis (km)
        e: eccentricity
        inc: inclination (radians)
        RAAN: right-ascension (radians)
        argPer: argument of perigee (radians)
        trueAnom: true anomaly (radians)
"""
function rv2coe(vec, mu) 
    rvec = view(vec, 1:3)
    vvec = view(vec, 4:6)

    r  =norm(rvec)
    v = norm(vvec)
    
    # ang. momentum
    hvec = cross(rvec, vvec)
    h = norm(hvec)
    
    I = SA[1, 0, 0]
    J = SA[0, 1, 0]
    K = SA[0, 0, 1]

    nvec = cross(K,hvec)
    n = norm(nvec)

    # eccentricity
    evec = ((v^2-mu/r)*rvec-dot(rvec,vvec)*vvec)/mu 
    e = norm(evec)

    #Semi-major and semi-latus(not returned for now)
    zeta = v^2/2-mu/r
    if e!=1.0
        a = -mu/(2*zeta)
        p = a*(1-e^2)
    else
        a=Inf
        p = h^2/mu
    end

    # Incl.
    i = acos(hvec[3]/h)

    # RAAN
    Ω = acos(nvec[1]/n)
    if nvec[2] < 0-eps()
        Ω = 2*pi-Ω
    end

    # Arg. per

    # I had some numerical issues with the acosarg being slightly (on order of 1e-10) higher than 1.0, causing a domain error. This is a band aid.
    acosarg = dot(nvec, evec) / (n*e)
    if acosarg > 1.0 || acosarg < -1.0
        clampval = clamp(acosarg, -1.0, 1.0)
        @warn("Argument of inverse cosine ($acosarg) out of domain for ω calculation. Clamping to $clampval. Answers may not be accurate.")
        acosarg = clampval
    end

    ω = acos(acosarg)
    if evec[3] < 0-eps()
        ω = 2*pi-ω
    end

    # T.A
    acosarg = dot(evec, rvec) / (e*r)
    if acosarg > 1.0 || acosarg < -1.0
        clampval = clamp(acosarg, -1.0, 1.0)
        @warn("Argument of inverse cosine ($acosarg) out of domain for θ calculation. Clamping to $clampval. Answers may not be accurate.")
        acosarg = clampval
    end
    θ = acos(acosarg)
    if dot(r,v)<0
        θ = 2*pi-θ
    end

    # special cases
    # i am lazy :) (pg 116 in vallado, do this at some point)
    
    return SA[a, e, i, Ω, ω, θ]
    
end