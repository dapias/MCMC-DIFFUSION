include("../src/MCMCDiffusion.jl")

using MCMCDiffusion, DynamicalBilliards, JLD
T = Float64

space = T(2.2)
r = T(1.0)

polygon_sides = n = 6
bth = billiard_polygon(polygon_sides, space/sqrt(3); setting = "periodic")
d = Disk(([T(0.),T(0.)]),r)
push!(bth, d)

D = 0.17
beta = -1.
N = 10
to = collect(10:5.25.)
means = zeros(to)
for i in 1:length(to)
    ts = to[i]
    rstar = ts/2.
    chaotic_results, acceptance = rMCMC(ts, N, beta, D)
    res = sum(exp.(beta*chaotic_results[chaotic_results .> rstar]))/N
    means[i] = res
    println(estimator)
end
                      
save("meanslorentzmcmc.jld", "means", means, "time", to)
