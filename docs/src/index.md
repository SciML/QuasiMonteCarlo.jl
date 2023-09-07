# QuasiMonteCarlo.jl: Quasi-Monte Carlo (QMC) Samples Made Easy

QuasiMonteCarlo.jl is a lightweight package for generating Quasi-Monte Carlo (QMC) samples
using various different methods.

## Installation

To install QuasiMonteCarlo.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("QuasiMonteCarlo")
```

## Get Started

### Basic API

```julia
using QuasiMonteCarlo, Distributions
lb = [0.1,-0.5]
ub = [1.0,20.0]
n = 5
d = 2

s = QuasiMonteCarlo.sample(n,lb,ub,GridSample())
s = QuasiMonteCarlo.sample(n,lb,ub,Uniform())
s = QuasiMonteCarlo.sample(n,lb,ub,SobolSample())
s = QuasiMonteCarlo.sample(n,lb,ub,LatinHypercubeSample())
s = QuasiMonteCarlo.sample(n,lb,ub,LatticeRuleSample())
s = QuasiMonteCarlo.sample(n,lb,ub,HaltonSample())
```

The output `s` is a matrix, so one can use things like `@uview` from
[UnsafeArrays.jl](https://github.com/oschulz/UnsafeArrays.jl) for a stack-allocated
view of the `i`th point:

```julia
using UnsafeArrays
@uview s[:,i]
```

### MC vs QMC

We illustrate the gain of QMC methods over plain Monte Carlo using the 5-dimensional example from Section 15.9 in the [book by A. Owen](https://artowen.su.domains/mc/qmcstuff.pdf).

```@example MCvsQMC; continued = true
f₁(𝐱) = prod(1 + √(12)/5*(xⱼ - 1/2) for xⱼ ∈ 𝐱)
μ_exact = 1 # = ∫ f₁(𝐱) d⁵𝐱
```

One can estimate the integral $\mu$ using plain Monte Carlo, or Quasi Monte Carlo or Randomized Quasi Monte Carlo. See the other section of this documentation for more information on the functions used in the example.

```@example MCvsQMC
using QuasiMonteCarlo, Random, Distributions
using Plots; default(fontfamily="Computer Modern")
Random.seed!(1234)
d = 5 # Dimension (= prime base for Faure net)
b = 2 # Base for Sobol net
m_max = 19
m_max_Faure = 8
N = b^m_max

# Generate sequences
seq_MC = QuasiMonteCarlo.sample(N, d, Uniform()) # Monte Carlo i.i.d Uniform sampling
seq_QMC_Sobol = QuasiMonteCarlo.sample(N, d, SobolSample()) # Sobol net
seq_RQMC_Sobol = QuasiMonteCarlo.sample(N, d, SobolSample(R = OwenScramble(base = b, pad = 32))) # Randomized version of Sobol net
seq_RQMC_Faure = QuasiMonteCarlo.sample(d^m_max_Faure, d, FaureSample(R = OwenScramble(base = d, pad = 32))) # Randomized version of Faure net

# Estimate the integral for different n with different estimator μ̂ₙ
μ_MC = [mean(f₁(x) for x in eachcol(seq_MC[:, 1:b^m])) for m in 1:m_max]
μ_QMC_Sobol = [mean(f₁(x) for x in eachcol(seq_QMC_Sobol[:, 1:b^m])) for m in 1:m_max]
μ_RQMC_Sobol = [mean(f₁(x) for x in eachcol(seq_RQMC_Sobol[:, 1:b^m])) for m in 1:m_max]
μ_RQMC_Faure = [mean(f₁(x) for x in eachcol(seq_RQMC_Faure[:, 1:d^m])) for m in 1:m_max_Faure]

# Plot the error |μ̂-μ| vs n
plot(b.^(1:m_max), abs.(μ_MC .- μ_exact), label="MC")
plot!(b.^(1:m_max), abs.(μ_QMC_Sobol .- μ_exact), label="QMC Sobol")
plot!(b.^(1:m_max), abs.(μ_RQMC_Sobol .- μ_exact), label="RQMC Sobol")
plot!(d .^(1:m_max_Faure), abs.(μ_RQMC_Faure .- μ_exact), label="RQMC Faure")
plot!(n -> n^(-1/2), b.^(1:m_max), c = :black, s = :dot, label = "n^(-1/2)")
plot!(n -> n^(-3/2), b.^(1:m_max), c = :black, s = :dash, label = "n^(-3/2)") 
# n^(-3/2) is the theoretical scaling for scrambled nets e.g. Theorem 17.5. in https://artowen.su.domains/mc/qmcstuff.pdf
xlims!(1, 1e6)
ylims!(1e-9, 1)
xaxis!(:log10)
yaxis!(:log10)
xlabel!("n", legend = :bottomleft)
ylabel!("|μ̂-μ|")
```

## Adding a new sampling method

Adding a new sampling method is a two-step process:

1. Add a new SamplingAlgorithm type.
2. Overload the sample function with the new type.

All sampling methods are expected to return a matrix with dimension `d` by `n`, where `d` is the dimension of the sample space and `n` is the number of samples.

**Example**

```julia
struct NewAmazingSamplingAlgorithm{OPTIONAL} <: SamplingAlgorithm end

function sample(n,lb,ub,::NewAmazingSamplingAlgorithm)
    if lb isa Number
        ...
        return x
    else
        ...
        return reduce(hcat, x)
    end
end
```

## Contributing

- Please refer to the
  [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
  for guidance on PRs, issues, and other matters relating to contributing to SciML.
- There are a few community forums:
  - the #diffeq-bridged channel in the [Julia Slack](https://julialang.org/slack/)
  - [JuliaDiffEq](https://gitter.im/JuliaDiffEq/Lobby) on Gitter
  - on the [Julia Discourse forums](https://discourse.julialang.org)
  - see also [SciML Community page](https://sciml.ai/community/)

## Reproducibility

```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>and using this machine and Julia version.</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```

```@example
using Pkg # hide
Pkg.status(;mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

```@raw html
You can also download the
<a href="
```

```@eval
using TOML
version = TOML.parse(read("../../Project.toml",String))["version"]
name = TOML.parse(read("../../Project.toml",String))["name"]
link = "https://github.com/SciML/"*name*".jl/tree/gh-pages/v"*version*"/assets/Manifest.toml"
```

```@raw html
">manifest</a> file and the
<a href="
```

```@eval
using TOML
version = TOML.parse(read("../../Project.toml",String))["version"]
name = TOML.parse(read("../../Project.toml",String))["name"]
link = "https://github.com/SciML/"*name*".jl/tree/gh-pages/v"*version*"/assets/Project.toml"
```

```@raw html
">project</a> file.
```
