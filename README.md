# HBQL

R code implementing the high-dimensional Bayesian Q-learning (HBQL) method proposed in the manuscript “Bayesian Variable Selection for Dynamic Treatment Regimes with High-Dimensional Covariates”.

This repository provides the core implementation of the proposed HBQL method and an illustrative simulation example.

## Files

- `gib.sig.R`: Function for sampling $\sigma^2$ in the Gibbs sampler.
- `theta.initial.R`: Function for initializing theta.
- `gib.theta.R`: Function for sampling theta in the Gibbs sampler.
- `gib.gamma.R`: Function for sampling gamma in the Gibbs sampler.
- `simulate_data.R`: Function for generating simulated datasets.
- `fun_select_gamma.R`: Function for selecting active covariates at stages 1 and 2 using HBQL.
- `Simulation_example.R`: Simulation example for the HBQL method. This script illustrates Setting 1 with `p1 = 302` and `p2 = 603`. Settings 2 and 3 can be implemented by modifying the corresponding model parameters specified in the script.
