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
    D = 0.17
    beta = -1.

    rstar = ts/2.
    chaotic_results, acceptance = MCMCDiffusion.rMCMC(ts, N, bth,  beta, D)
    res = sum(exp.(beta*chaotic_results[chaotic_results .> rstar]))/N
end

to = collect(10.:2.5:22.5)
means = pmap(estimation, to)
save("meanslorentzmcmc.jld", "means", means, "time", to)
