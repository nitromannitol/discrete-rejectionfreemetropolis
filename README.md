# discrete-rejectionfreemetropolis

Julia code for a discrete rejection-free Metropolis algorithm, compared with the usual
Metropolis algorithm. The state is a configuration of `k` particles on a periodic `n × n`
grid, with the constraint that particles are at periodic Manhattan distance greater than a
given distance from each other. A rejection-free move starts at a particle and pushes along
`right` or `up` with a displacement parameter `λ`, following a chain of collisions
(`runCollisionChain` in [`functions.jl`](functions.jl)). The code builds the transition
matrices of the rejection-free chain and of the usual Metropolis chain for small `n` and `k`
(`generateTransitionMatrix`, `generateTransitionMatrixMetropolis`), lists all possible
configurations, and compares the second largest eigenvalues of the two chains. The figures in
[`figures/`](figures/) are produced by [`script.jl`](script.jl).

## Running

The code was written in March 2016 for a Julia release before 1.0 and does not run on
current Julia: `script.jl` uses `eig`, which was removed in Julia 1.0, and `norm` is called
without loading `LinearAlgebra`. The only package the code loads is `PyPlot` (matplotlib
through PyCall).

With a pre-1.0 Julia and `PyPlot` installed, run from the repository directory:

```julia
include("script.jl")
```

The script runs a sample chain, lists all configurations for `n = 2, k = 2` and `n = 3, k =
3`, shows configurations and moves that give invalid chains for `n = 4` and `n = 5` when `λ`
is fixed, prints the second largest eigenvalues of the two chains for `n = 2, …, 5`, and
plots the comparison of the second largest eigenvalue against `k` from values stored in the
script. The plots are saved as EPS files in `figures/`.

## Repository layout

| Path | Contents |
|---|---|
| [`functions.jl`](functions.jl) | configurations, collision chains, transition matrices, plotting |
| [`script.jl`](script.jl) | the runs that produce the figures |
| [`figures/`](figures/) | the figures as EPS files, with PDF conversions |

## Authors and citation

The code is by **Ahmed Bou-Rabee**.
