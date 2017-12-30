export Proposal, ShiftProposal, LocalProposal, shift_proposal, local_proposal

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

function local_proposal(init::InitialCondition, n::Int, bt::Vector{<:Obstacle{T}}, sigma::T) where {T <: AbstractFloat}

    delta_theta = rand()*(2pi)
    rprime = abs(randn()*sigma)

    s_new = cos(delta_theta)*rprime + init.s
    s_new = mod(s_new, 1.0)
    sinphi_new = sin(delta_theta)*rprime + init.sinphi

    if abs(sinphi_new) > 1.0
        return init
    end

    p = particle_from_coordinates([s_new, sinphi_new], n, bt)
    side = ceil(Int, s_new*n)

    return InitialCondition(p, s_new, sinphi_new, side)
end

"""
Function that evolves a particle until it reaches a periodic wall. It returns the index associated with the periodic image of the side that the particle reaches (where the motion would continue)
"""
function evolution_to_periodicwall!(p::Particle{T}, bt::Vector{<:Obstacle{T}}, n::Int) where {T<: AbstractFloat}
    i = 1
    while true
        tmin::T, i::Int = next_collision(p, bt)
        tmin = DynamicalBilliards.relocate!(p, bt[i], tmin)
        resolvecollision!(p, bt[i])
        if typeof(bt[i]) <: PeriodicWall
            break
        end
    end
    if n == 6
        index = i >= 4 ? mod(i + div(n,2), n) : i + div(n,2)
    elseif n==4
        index = i >= 3 ? mod(i + div(n,2), n) : i + div(n,2)
    end
        
end

"""
Function that is passed to the type `ShiftProposal`
"""
function shift_proposal(particle::Particle{T}, sides::Int, bt::Vector{<:Obstacle{T}}, tshift::T, index::Int) where {T <: AbstractFloat}
    p = copy(particle)
    if tshift < 0.0  ##Invert the direction of the velocity
        #and put the particle on its periodic image
        p.vel = -p.vel
        p.pos += bt[index].normal 
    end
    
    t = abs(tshift)
    count = zero(T)
    
    while count < t
        tmin::T, i::Int = next_collision(p, bt)
        # set counter
        if count +  DynamicalBilliards.increment_counter(t, tmin) > t
            break
        end
        #####
        tmin = DynamicalBilliards.relocate!(p, bt[i], tmin)
        resolvecollision!(p, bt[i])
        count += DynamicalBilliards.increment_counter(t, tmin)
    end#time loop
    ##Here the evolution ends and then we propagate until the next collision
    new_index = evolution_to_periodicwall!(p, bt, sides)
    new_particle = Particle(p.pos, p.vel, SVector{2,T}([0.,0.]))

    if tshift < 0.0
        new_particle.vel = -new_particle.vel
        new_particle.pos += bt[new_index].normal
        if sides == 6
            new_index = new_index >= 4 ? mod(new_index + div(sides,2), sides) : new_index + div(sides,2)
        elseif sides == 4
            new_index = new_index >= 3 ? mod(new_index + div(sides,2), sides) : new_index + div(sides,2)
        end
    end
    s, sinphi = coordinates_from_particle(new_particle, sides, bt, new_index)
    
    return InitialCondition(new_particle, s, sinphi, new_index)
end




