export rMCMC


function tstar(to::S, beta::Float64,  D::Float64, y::T; ac = 0.5) where {T <: AbstractFloat, S<: Integer}
    tst = to - abs((ac-1)/(2*D*beta) * 1/(1 - y^2))
    if tst < to*9/10
        return T(to*9/10)
    end
    tst
end

function sigma(x::T, to::S, y::T, beta::Float64, D::Float64; a = 4.) where {T <: AbstractFloat, S<: Integer}
    lambda = log(a)
    tst = tstar(to, beta, D, y)
    exp(-lambda*tst)
end

function t_shift(to::S, beta::Float64, D::Float64, y::T) where {T<:AbstractFloat, S<: Integer}
    tst = tstar(to, beta, D, y)
    to - tst
end

function ratio_local(x1::T, x2::T, sigma1::T, sigma2::T) where {T<: AbstractFloat} 
    d  = abs(x1 - x2)
    expo = -d^2/2*(1./sigma2^2. - 1./sigma1^2.)
    return exp(expo)*sigma1/sigma2
end

function ratio_shift(tshift::T, tmean1::T, tmean2::T, sigma_t::Float64, to::S) where {T<:AbstractFloat, S<:Integer}
    
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

#    f1 = phi(xi1)/sigma_t
#    f2 = phi(xi2)/sigma_t

    return f2/f1
end


function proposal(s1::ShiftProposal, s2::ShiftProposal, l::LocalProposal)
    r = rand()
    if r < 0.5
        if r < 0.25
            return s1
        else
            return s2
        end
    else
        return l
    end
end

function proposal(tshift, sigma, a)
    shift1 = ShiftProposal(tshift, x -> shift_proposal(x, tshift, a))
    shift2 = ShiftProposal(-tshift, x -> shift_proposal(x, -tshift, a))
    local_prop = LocalProposal(sigma, x -> local_proposal(x, sigma) )
    
    proposal(shift1, shift2, local_prop)
end

function proposal(forw_prop::ShiftProposal, sigma, a)
    tshift = -forw_prop.parameter
    ShiftProposal(tshift,
                  x -> shift_proposal(x, tshift, a))
end

function proposal(forw_prop::LocalProposal, sigma, a)
    LocalProposal(sigma, x -> local_proposal(x, sigma) )
end

function ratio(x1::T, x2::T, sigma1::T, sigma2::T, tshift::T, tmean1::T, tmean2::T, sigma_t::Float64, to::S, forward::ShiftProposal, backward::ShiftProposal) where {T <: AbstractFloat, S <: Integer}

    ratio_shift(tshift, tmean1, tmean2, sigma_t, to)

end

function ratio(x1::T, x2::T, sigma1::T, sigma2::T, tshift::T, tmean1::T, tmean2::T, sigma_t::Float64, to::S, forward::LocalProposal, backward::LocalProposal)  where {T <: AbstractFloat, S <: Integer}

    ratio_local(x1, x2, sigma1, sigma2)
end


function rMCMC(to::S, N::Int,beta::Float64, D::Float64; T = Float64, sigma_t = 1.0, a_parameter = 4.) where {S<: Integer}
    birk_coord = zeros(T, N,2)
    ###initialize
    init  = init_prime = T(rand())
    birk_coord[1, 1] = init
    y_prime = y = abs(evolution(init, to, a_parameter) - 1/2.)/sqrt(2*D*to)
    birk_coord[1, 2] = y
    acceptance = 0
    ######
    for i in 2:N
        tshift_mean = t_shift(to, beta, D, y)
        tshift = T(rand(TruncatedNormal(tshift_mean, sigma_t, 0., to)))
        sigma_local =  sigma(init, to, y ,beta, D, a = a_parameter)
        forw_prop = proposal(tshift, sigma_local, a_parameter)
        
        init_prime = forw_prop.f(init)
        y_prime = abs(evolution(init_prime, to, a_parameter) - 1/2.)/sqrt(2*D*to)
        tshift_meanprime = t_shift(to, beta, D, y_prime)
        sigma_prime = sigma(init_prime, to, y_prime, beta, D, a = a_parameter)
        back_prop = proposal(forw_prop, sigma_prime, a_parameter)
        
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
        
        birk_coord[i, :] = [init, y*sqrt(2*D*to)]
    end
    
    
    birk_coord, acceptance/N
end

