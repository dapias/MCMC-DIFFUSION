__precompile__()

module MCMCDiffusion

using DynamicalBilliards, Distributions

include("Boxmap.jl")
include("Proposals.jl")
include("MCMC.jl")

end
