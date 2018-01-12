export randominitialcondition, InitialCondition


"""
Type that contains the particle, its Birkhoff coordinates and an index representing the side of the polygon where it is located
"""
struct InitialCondition{T<:AbstractFloat}
    particle::Particle{T}
    phi::T
end
"""
Returns the tuple `(p, s, phi)` being `p` of type `Particle` located randomly in the boundary of the unit cell with Birkhoff coordinates `s`, `phi`. If sides is true it also returns the side of the cell where the particle is located
"""
function randominitialcondition(bt::Vector{<:Obstacle{T}}) where {T<:AbstractFloat}

    xmin::T, ymin::T, xmax::T, ymax::T = DynamicalBilliards.cellsize(bt)
    f = T(rand())
    while f == 0 || f==1/4 || f==1/2 || f == 3/4
        f = T(rand())
    end
    phi0 = f*2*pi - pi   ##phi0 between -pi and pi
    xp = T(rand())*(xmax-xmin) + xmin
    yp = T(rand())*(ymax-ymin) + ymin
    p = Particle([xp, yp, phi0])
    
    dist = DynamicalBilliards.distance_init(p, bt)
    while dist <= sqrt(eps(T))

        xp = T(rand())*(xmax-xmin) + xmin
        yp = T(rand())*(ymax-ymin) + ymin
        p.pos = SVector{2,T}(xp, yp)
        dist = DynamicalBilliards.distance_init(p, bt)
    end

    init = InitialCondition(p, phi0)

    return init
end
