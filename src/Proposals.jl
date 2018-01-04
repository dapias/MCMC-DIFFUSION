export Proposal, ShiftProposal, LocalProposal, shift_proposal, local_proposal, billiard_evolution!
import Base.copy 

function copy(p::Particle{T}) where {T<:AbstractFloat}
    pos = copy(p.pos)
    vel = copy(p.vel)
    cc = copy(p.current_cell)

    Particle(pos,vel, cc)
end



abstract type Proposal end

mutable struct ShiftProposal{T} <: Proposal
    parameter::Union{Function, T}    #t_shift
    f::Function
end

mutable struct LocalProposal{T} <: Proposal
    parameter::Union{Function, T}    #Standard_Deviation in the case of the Gaussian
    f::Function
end

##For Box map

function shift_proposal(x::T, tshift::T, a::Float64) where {T<:AbstractFloat}
    tk = 1
    if tshift > 0.
        tk = Int(ceil(tshift))
    else
        tk = Int(floor(tshift))
    end
    xprime = evolution(x, tk, a)  

    xprime = mod(xprime, 1.)  #Translation of the box
    return xprime
end


function local_proposal(x::T, sigma::T) where {T<: AbstractFloat}
   
    xprime = T(randn())*sigma + x

    if xprime < 0.0 || xprime > 1.0  ##Prohibited point
        return x
    end
    
    xprime
end



##For Billiards
function billiard_evolution!(p::Particle{T}, bt::Vector{<:Obstacle{T}}, t::T) where {T<:AbstractFloat}

    count = zero(T)

    while count < t
        tmin::T, i::Int = next_collision(p, bt)
        # set counter
        if count +  DynamicalBilliards.increment_counter(t, tmin) > t
            break
        end
        #####
        tmin = relocate!(p, bt[i], tmin)
        resolvecollision!(p, bt[i])
        count += DynamicalBilliards.increment_counter(t, tmin)
    end#time loop

    tmin = t - count 
    propagate!(p, tmin)
    
    return nothing
end

function local_proposal(init::InitialCondition, bt::Vector{<:Obstacle{T}}, sigma::T) where {T <: AbstractFloat}

    rprime = abs(randn()*sigma)
    rand_theta = T(rand())*2*pi

    d1 = abs(rprime*cos(rand_theta))
    d2 = abs(rprime*sin(rand_theta))

    xmin::T, ymin::T, xmax::T, ymax::T = DynamicalBilliards.cellsize(bt)

    pnew = copy(init.particle)

    theta = T(rand())*(2pi)
    
    xnew = d1*cos(theta) + pnew.pos[1]
    ynew = d1*sin(theta) + pnew.pos[2]
    rdir = rand([-1,1])
    phinew = mod(init.phi + rdir*d2, 2pi)

    if (xmin <= xnew <= xmax) && (ymin <= ynew <= ymax)
        pnew = Particle([xnew, ynew, phinew])
        dist = DynamicalBilliards.distance_init(pnew, bt)
        new_init = InitialCondition(pnew, phinew)
        
        if dist <= sqrt(eps(T))
            return init
        else
            return new_init
        end
    else
        return init
    end
end

"""
Function that is passed to the type `ShiftProposal`
"""
function shift_proposal(init::InitialCondition, bt::Vector{<:Obstacle{T}}, tshift::T) where {T <: AbstractFloat}

    p = copy(init.particle)
    if tshift < 0.0  ##Invert the direction of the velocity
        p.vel = -p.vel
    end
    
    t = abs(tshift)
    billiard_evolution!(p, bt, t)

    if tshift < 0.0
        p.vel = -p.vel
    end
    

    phinew = atan2(p.vel[2], p.vel[1])
    if phinew < 0.0
        phinew += 2pi
    end

    p = Particle([p.pos[1],p.pos[2], phinew])
    init = InitialCondition(p, phinew)

    return init
end




