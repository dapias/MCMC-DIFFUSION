include("../src/MCMCDiffusion.jl")

using JLD, MCMCDiffusion, DynamicalBilliards, LsqFit

function rdirect(to, n)
    distances = zeros(n);
    for i in 1:n
        p = randominitialcondition(bt).particle
        distances[i] = MCMCDiffusion.distance(p, bth, to)
    end
    distances
end

function costdirect(to, nsimulations, mean_value; prec = 10.^-1., D = 0.17)
    ystar = to/2.
    N = Int.(logspace(3,6, 4))
    i =1
    estimator_out = zeros(length(N))
    estimator_in = zeros(nsimulations)
    n = N[i]
    for j in 1:nsimulations
        res = rdirect(to, n)
        estimator_in[j] =  length(res[res .>= ystar])/n
    end
    estimator_out[1] = std(estimator_in, mean = mean_value)
    for i in 2:length(N)
        n = N[i]
        for j in 1:nsimulations
            res = rdirect(to, n) 
            estimator_in[j] =  (length(res[res .>= ystar])/n)
        end
        estimator_out[i] = std(estimator_in, mean = mean_value)  
    end
   
    
    model(x, p) = p[1] + p[2]*x
    xdata = log.(N)
    ydata =  log.(estimator_out/mean_value)
    p0 = [1.,-1.e-3]
    fit = curve_fit(model, xdata, ydata, p0)
    estimator =  exp((log(prec) - fit.param[1])/fit.param[2])
    
end

function costdirect(to, nsimulations, ngroups, mean_value)
    estimator = zeros(ngroups)
    for i in 1:ngroups
        estimator[i] = costdirect(to, nsimulations, mean_value)
    end
    estimator
end

function estimatedirect(tmin, tmax, nsimulations, ngroups, means)
    t = collect(tmin:5.:tmax)
    i = 1
    mean_value = means[i]
    est = costdirect(t[i], nsimulations, ngroups, mean_value)
    result = est

    for i in 2:length(t)
        mean_value = means[i]
        est = costdirect(t[i], nsimulations, ngroups, mean_value)
        result = hcat(result, est)
    end

    filename = "directcostlorentz"*randstring(2)*".jld"
    save(filename, "time", t, "estimator", result)
    println("Filename $filename generated")
end


T = Float64

space = T(2.2) # Space between two adjacent disks
r = T(1.0) # Radius of the disk
polygon_sides = 6
bth = bt =  billiard_polygon(polygon_sides, space/sqrt(3); setting = "periodic")
d = Disk(([T(0.),T(0.)]), r)
push!(bth, d);

a = load("meanslorentzdirect.jld")
means = a["means"]

# estimatedirect(10., 30., 10, 1, means)
