include("../src/MCMCDiffusion.jl")

using MCMCDiffusion, DynamicalBilliards, JLD
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
        p = randominitialcondition(bt).particle
        distances[i] = MCMCDiffusion.distance(p,bth,to)
    end
    distances
end

N = 10
to = collect(10:5.25.)
means = zeros(to)
for i in 1:length(to)
    ts = to[i]
    rstar = ts/2.
    res = rdirect(ts, N)
    estimator = length(res[res .>= rstar])/N
    means[i] = estimator
    println(estimator)
end
                      
save("meanslorentzdirect.jld", "means", means, "time", to)
