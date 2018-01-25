addprocs(5)

@everywhere include("../src/MCMCDiffusion.jl")
@everywhere using MCMCDiffusion, DynamicalBilliards, JLD


@everywhere function estimation(ts)
    N = 10^7
    
    T = Float64
    space = T(2.2)
    r = T(1.0)
    polygon_sides = n = 6
    bth = billiard_polygon(polygon_sides, space/sqrt(3); setting = "periodic")
    d = Disk(([T(0.),T(0.)]),r)
    push!(bth, d)

    function rdirect(to, n)
        distances = zeros(n)
        for i in 1:n
            p = MCMCDiffusion.randominitialcondition(bth).particle
            distances[i] = MCMCDiffusion.distance(p,bth,to)
        end
        distances
    end

    
    rstar = ts/2.
    res = rdirect(ts, N)
    estimator = length(res[res .>= rstar])/N
end

to = collect(10.:2.5:22.5)
means = pmap(estimation, to)
save("meanslorentzdirect.jld", "means", means, "time", to)
