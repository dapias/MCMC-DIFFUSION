export rMCMC, distance


function tstar(beta::Float64, to::T, D::Float64, y::T; a = 0.5) where {T <: AbstractFloat}
    tstar = to - abs((a-1.)/(2*D*beta) * 1/(1-y^2))
    if tstar < to*9/10
        return to*9/10
    end
    tstar
end

function sigma(beta::Float64, t::T, D::Float64, y::T; Delta = T(1.0),
               l_star = 2.10)  where {T <: AbstractFloat}
    t_star =  tstar(beta, t, D, y)
    Delta*exp(-l_star*t_star)
end

function t_shift(beta::Float64, to::T, D::Float64, y::T) where {T<: AbstractFloat}
    tst = tstar(beta, to, D, y)
    to - tst
end

function ratio_local(init1::InitialCondition{T}, init2::InitialCondition{T}, sigma1::T, sigma2::T) where {T <: AbstractFloat}
    d  = distance(init1, init2)
    expo = -d^2/2*(1./sigma2^2. - 1./sigma1^2.)
    return exp(expo)*sigma1/sigma2
end

##Using the product metric space 
function distance(i1::InitialCondition{T}, i2::InitialCondition{T}) where {T<:AbstractFloat}
    x1 = [i1.particle.pos[1], i1.particle.pos[2]]
    x2 = [i2.particle.pos[1], i2.particle.pos[2]]

    d1 = norm(x1 -x2)

    phi1 = i1.phi
    phi2 = i2.phi
    deltaphi = abs(phi2 -phi1)
    if deltaphi > pi
        deltaphi = 2pi - abs(deltaphi)
    end
    new_vector = [d1, deltaphi]

    norm(new_vector)
end

function proposal(tshift::T, sigma::T, bt::Vector{<:Obstacle{T}}) where {T <: AbstractFloat}
    shift1 = ShiftProposal(tshift, x::InitialCondition -> shift_proposal(x, bt, tshift))
    shift2 = ShiftProposal(-tshift, x::InitialCondition -> shift_proposal(x, bt, -tshift))
    local_prop = LocalProposal(sigma, x::InitialCondition -> local_proposal(x, bt, sigma) )

    proposal(shift1, shift2, local_prop)
end
    
function ratio(x1::InitialCondition{T}, x2::InitialCondition{T}, sigma1::T, sigma2::T, tshift::T, tmean1::T, tmean2::T, sigma_t::Float64, to::T, forward::ShiftProposal, backward::ShiftProposal) where {T <: AbstractFloat}

    ratio_shift(tshift, tmean1, tmean2, sigma_t, to)

end

function ratio_shift(tshift::T, tmean1::T, tmean2::T, sigma_t::Float64, to::T) where {T<:AbstractFloat}
    
    xi1 = (tshift - tmean1)/sigma_t
    xi2 = (tshift - tmean2)/sigma_t

    alpha1 = (0 - tmean1)/sigma_t
    beta1 = (to - tmean1)/sigma_t
   
    alpha2 = ( 0 - tmean2)/sigma_t
    beta2 = (to - tmean2)/sigma_t

    phi(x) = 1/sqrt(2*pi)*exp(-1/2*x^2)
    psi(x) = 1/2*(1 + erf(x/sqrt(2)))

    f1 = phi(xi1)/(sigma_t*(psi(beta1) - psi(alpha1)))
    f2 = phi(xi2)/(sigma_t*(psi(beta2) - psi(alpha2)))

    return f2/f1
end



function ratio(x1::InitialCondition{T}, x2::InitialCondition{T}, sigma1::T, sigma2::T, tshift::T, tmean1::T, tmean2::T, sigma_t::Float64, to::T, forward::LocalProposal, backward::LocalProposal)  where {T <: AbstractFloat}

    ratio_local(x1, x2, sigma1, sigma2)
end


function proposal(forw_prop::ShiftProposal, sigma::T, bt::Vector{<:Obstacle{T}}) where {T <: AbstractFloat}
    tshift = -forw_prop.parameter
    ShiftProposal(tshift, x::InitialCondition -> shift_proposal(x, bt, tshift))
end

function proposal(forw_prop::LocalProposal, sigma, bt::Vector{<:Obstacle{T}}) where {T <: AbstractFloat}
    LocalProposal(sigma, x::InitialCondition -> local_proposal(x, bt, sigma) )
end

function distance(particle::Particle{T}, bt::Vector{<:Obstacle{T}}, t::T) where {T<:AbstractFloat}
    p = copy(particle)
    rpos = SVector{2,T}[]
    push!(rpos, p.pos)
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
    push!(rpos, p.pos + p.current_cell)
    d = norm(rpos[2]-rpos[1])
end

function rMCMC(to::T, N::Int, bt::Vector{<:Obstacle{T}}, beta::Float64, D::Float64; sigma_t = 2.0) where {T<: AbstractFloat}

    birk_coord = zeros(T, N, 1)
    ###initialize
    init  = init_prime = randominitialcondition(bt)
    obs = x::Particle -> distance(x, bt, to)/sqrt(2*D*to) ##y-observable
    y_prime = y = obs(init.particle)  ##Actually this is y
    birk_coord[1] = y*sqrt(2*D*to)
    acceptance = 0
    ####################
    for i in 2:N
        tshift_mean = t_shift(beta, to, D, y)   
        tshift = T(rand(Normal(tshift_mean, sigma_t)))
        sigma_local =  sigma(beta, to, D, y)
        forw_prop = proposal(tshift, sigma_local, bt)

        init_prime = forw_prop.f(init)
        y_prime = obs(init_prime.particle)
        tshift_meanprime = t_shift(beta, to, D, y_prime)
        sigma_prime =  sigma(beta, to, D, y_prime)
        back_prop = proposal(forw_prop, sigma_prime, bt)

        rat = ratio(init, init_prime, sigma_local, sigma_prime, tshift, tshift_mean, tshift_meanprime, sigma_t, to, forw_prop, back_prop)
        ac = rat*exp(-beta*(y_prime -y)*sqrt(2*D*to))

        if ac >= 1.0
            init = init_prime
            y = y_prime
            acceptance += 1
        else
            rad = rand()
            if rad < ac
                init = init_prime
                y = y_prime
                acceptance += 1
            end
        end
        
        birk_coord[i] = y*sqrt(2*D*to)
    end
    
    birk_coord, acceptance/(N-1)
end




