##The functions that defines the dynamics are defined. Also, the functions for finding the preimages.

function boxmap(x::T, a::Float64) where {T <: AbstractFloat}
    if 0. <= x <= 1/2.
        return  a*x
    elseif 1/2. < x <= 1.
        return a*x+ 1.-a
    end
end

function Ma(x::T, a::Float64) where {T<: AbstractFloat}
   if 0. <= x <= 1.
        return boxmap(x, a)
    elseif x > 1.
        return Ma(x-1.,a) + 1.
    elseif x < 0.
        return -Ma(-x, a)
    end
end

function inverseboxmap(x::T, a::Float64) where {T <: AbstractFloat}  
    if x ==0.0
        r = rand()
        if r < 1/3.
            return x/a
        elseif 1/3. < r < 2/3.
            return (x + a - 1)/a
        else
            return (x - a + 1)/a
        end
        
    elseif x == 1.
        r = rand()
        if r < 1/3.
            return x/a
        elseif 1/3. < r < 2/3.
            return (x + a - 1)/a
        else
            return (x + a -2.)/a + 1
        end

    else
        r = rand()
        if r < 1/a
            return x/a
        elseif 1/a < r < 1/2.
            return (x + a - 1)/a
        elseif 1/2. < r < (1 - 1/a)
            return (x+ a -2.)/a + 1
        else
            return (x - a + 1)/a
        end
    end
end

function inverseMa(x::T, a::Float64) where {T <: AbstractFloat}

    if  0.<= x <= 1.
        return inverseboxmap(x, a)
    end


    if x > 1.
        modu = mod(x,1.)
        divi = div(x,1.)
        if modu == 0.0
            y = inverseboxmap(T(1.), a)
            return y + divi -1
        else
            y = inverseboxmap(modu, a)
            return y + divi
        end
        
    elseif x < 0.
        x = -x
        return -inverseMa(x, a)
    end    
end



function evolution(x::T, to::S, a::Float64) where {T<: AbstractFloat, S <: Integer}
    if to == 0
        return x
    elseif to > 0
        for k in 1:to
            x = Ma(x, a)
        end
    else 
        for k in 1:abs(to)
            x = inverseMa(x, a)
        end
    end
    x
end





    
