# Metropolis - Hastings for diffusive deterministic dynamical systems

We provide the code that supports the results reported in the manuscript [Monte Carlo sampling in diffusive dynamical systems](https://arxiv.org/abs/1804.06698). We implement here the Metropolis-Hastings algorithm for the estimation of a weighted distribution of the displacement.

In the folder `Examples` we illustrate its basic use for two dynamical systems: Lorentz gas and Box map. For its use in other deterministic dynamical systems the dynamics should be introduced and the function `rMCMC`. The steps to do it are always the same, namely: 

1) The dynamical system is coded, taking care of defining it for both positive and negative times. You might use the package [DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl) .

2) The main parameters, i.e. the mean Lyapunov exponent and the diffusion coefficient are estimated or passed if they are previously known.

3) The two proposals are coded based on the dynamics. It would depend on the nature of the dynamical system (e.g. map or flow) and its dimension.

4) The observable $r_{t_o}$ is defined and used inside `rMCMC`. In our examples, it is the `distance` for the Lorentz gas and `abs(evolution - 1/2)` for the box map. 



## Installation

From within Julia do

```
julia> Pkg.clone(https://github.com/dapias/MCMC-DIFFUSION.git")
```

 
