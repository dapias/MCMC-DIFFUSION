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





