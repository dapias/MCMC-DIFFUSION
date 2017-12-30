__precompile__()

module MCMCDiffusion

using DynamicalBilliards, Distributions

include("boxmap.jl")
include("proposals.jl")
include("MCMC.jl")

end
