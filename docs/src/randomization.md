# [Randomization methods](@id Randomization)

Most of the methods presented in [Sampler](@ref Samplers) are deterministic, i.e. `X = sample(n, d, ::DeterministicSamplingAlgorithm)` will always produce the same sequence $X = (X_1, \dots, X_n)$.

The main issue with deterministic Quasi Monte Carlo sampling is that it does not allow easy error estimation as opposed to plain Monte Carlo, where the variance can be estimated.

A Randomized Quasi Monte Carlo method must respect the two following criteria:

 1. Have $X_i\sim \mathbb{U}([0,1]^d)$ for each $i\in \{1,\cdots, n\}$.
 2. Preserve the QMC properties, i.e. the randomized $X$ still has low discrepancy.

This randomized version is unbiased and can be used to obtain a confidence interval or to do sensitivity analysis.

A good reference is the [book by A. Owen](https://artowen.su.domains/mc/qmcstuff.pdf), especially the Chapters 15, 16 and 17.

## API for randomization

```julia
abstract type RandomizationMethod end
```

There are two ways to obtain a randomized sequence:

- Either directly use `QuasiMonteCarlo.sample(n, d, DeterministicSamplingAlgorithm(R = SomeRandomizationMethod()))` or `sample(n, lb, up, DeterministicSamplingAlgorithm(R = RandomizationMethod()))`.
- Or, given $n$ points $d$-dimensional points, all in $[0,1]^d$ one can do `randomize(X, SomeRandomizationMethod())` where $X$ is a $d\times n$-matrix.

```@docs
randomize
```

The default method of `DeterministicSamplingAlgorithm` is `NoRand`

```@docs
NoRand
```

To obtain multiple independent randomizations of a sequence, i.e. Design Matrices, look at the [Design Matrices section](@ref DesignMatrices).

!!! note

    In most other QMC packages, randomization is performed "online" as the points are samples. Here, randomization is performed after the deterministic sequence is generated. Both methods are useful in different contexts. The former is generally faster to produce one randomized sequence, while the latter is faster to produce independent realizations of the sequence.
    
    **PRs are welcomed** to add "online" version of the sequence! See [this comment for inspiration](https://github.com/SciML/QuasiMonteCarlo.jl/pull/57#issuecomment-1326662016).
    
    Another way to view the two approaches is: given a computational budget of $N$ points, one can
    
     1. Put all of it into a sequence of size, $N$, thus having the best estimator $\hat{\mu}_N$. The price to pay is that this estimation is not associated with a variance estimation.
     2. Divide your computational budget into $N = n\times M$ to get $M$ independent estimator $\hat{\mu}_n$. From there, one can compute the empirical variance of the estimator.

## Scrambling methods

```julia
abstract type ScrambleMethod <: RandomizationMethod end
```

```@docs
ScrambleMethod
```

`ScramblingMethods(b, pad, rng)` are well suited for $(t,m,d)$-nets in base $b$. `b` is the base used to scramble and `pad` the number of bits in `b`-ary decomposition, i.e. $y \simeq \sum_{k=1}^{\texttt{pad}} y_k/\texttt{b}^k$.

The `pad` is generally chosen as $\gtrsim \log_b(n)$.

!!! warning

    In principle, the base `b` used for scrambling methods `ScramblingMethods(b, pad, rng)` can be an arbitrary integer.
    However, to preserve good Quasi Monte Carlo properties, it must match the base of the sequence to scramble.
    For example, (deterministic) Sobol sequences are base $b=2$, $(t,m,d)$ sequences while (deterministic) Faure sequences are $(t,m,d)$ sequences in prime base i.e. $b$ is an arbitrary prime number.

The implemented `ScramblingMethods` are

- `DigitalShift` the simplest and fastest method. For a point $x\in [0,1]^d$ it does $y_k = (x_k + U_k) \mod b$ where $U_k \sim \mathbb{U}(\{0, \cdots, b-1\})$

```@docs
DigitalShift
```

- `MatousekScramble` a.k.a. Linear Matrix Scramble is what people use in practice. Indeed, the observed performances are similar to `OwenScramble` for a lesser numerical cost.

```@docs
MatousekScramble
```

- `OwenScramble` a.k.a. Nested Uniform Scramble is the most understood theoretically, but is more costly to operate.

```@docs
OwenScramble
```

## Other methods

`Shift(rng)` a.k.a. Cranley-Patterson Rotation. It is by far the fastest method; it is used in `LatticeRuleScramble` for example.

```@docs
Shift
```

## Example

Randomization of a Faure sequence with various methods.

### Generation

```@example 1
using QuasiMonteCarlo, Random
Random.seed!(1234)
m = 4
d = 3
b = QuasiMonteCarlo.nextprime(d)
N = b^m # Number of points
pad = m # Can also choose something as `2m` to get "better" randomization

# Unrandomized low discrepancy sequence
x_faure = QuasiMonteCarlo.sample(N, d, FaureSample())

# Randomized version
x_uniform = rand(d, N) # plain i.i.d. uniform
x_shift = randomize(x_faure, Shift())
x_nus = randomize(x_faure, OwenScramble(base = b, pad = pad)) # equivalent to sample(N, d, FaureSample(R = OwenScramble(base = b, pad = pad)))
x_lms = randomize(x_faure, MatousekScramble(base = b, pad = pad))
x_digital_shift = randomize(x_faure, DigitalShift(base = b, pad = pad))
```

### Visualization of different methods

Plot the resulting sequences along dimensions `1` and `3`.

```@example 1
using Plots
# Setting I like for plotting
default(fontfamily = "Computer Modern",
    linewidth = 1,
    label = nothing,
    grid = true,
    framestyle = :default)

d1 = 1
d2 = 3 # you can try every combination of dimension (d1, d2)
sequences = [x_uniform, x_faure, x_shift, x_digital_shift, x_lms, x_nus]
names = [
    "Uniform",
    "Faure (deterministic)",
    "Shift",
    "Digital Shift",
    "Matousek Scramble",
    "Owen Scramble",
]
p = [plot(thickness_scaling = 1.5, aspect_ratio = :equal) for i in sequences]
for (i, x) in enumerate(sequences)
    scatter!(p[i], x[d1, :], x[d2, :], ms = 1.5, c = 1, grid = false)
    title!(names[i])
    xlims!(p[i], (0, 1))
    ylims!(p[i], (0, 1))
    yticks!(p[i], [0, 1])
    xticks!(p[i], [0, 1])
    hline!(p[i], range(0, 1, step = 1 / 4), c = :gray, alpha = 0.2)
    vline!(p[i], range(0, 1, step = 1 / 4), c = :gray, alpha = 0.2)
    hline!(p[i], range(0, 1, step = 1 / 2), c = :gray, alpha = 0.8)
    vline!(p[i], range(0, 1, step = 1 / 2), c = :gray, alpha = 0.8)
end
plot(p..., size = (800, 600))
```

### $(t,m,d)$-net visualization

Faure nets and their scrambled versions are digital $(t,m,d)$-net, which means they have strong equipartition properties.
On the following plot, we can (visually) verify that with Nested Uniform Scrambling, it also works with Linear Matrix Scrambling and Digital Shift.
You must see one point per rectangle of volume $1/b^m$. Points on the "left" border of rectangles are included, while those on the "right" are excluded. See [Chapter 15.7](https://artowen.su.domains/mc/qmcstuff.pdf) and Figure 15.10 for more details.

```@example 1
d1 = 1
d2 = 3 # you can try every combination of dimension (d1, d2)
x = x_nus # Owen Scramble, you can try x_lms and x_digital_shift
p = [plot(thickness_scaling = 1.5, aspect_ratio = :equal) for i in 0:m]
for i in 0:m
    j = m - i
    xᵢ = range(0, 1, step = 1 / b^(i))
    xⱼ = range(0, 1, step = 1 / b^(j))
    scatter!(p[i + 1], x[d1, :], x[d2, :], ms = 1.5, c = 1, grid = false)
    xlims!(p[i + 1], (0, 1.01))
    ylims!(p[i + 1], (0, 1.01))
    yticks!(p[i + 1], [0, 1])
    xticks!(p[i + 1], [0, 1])
    hline!(p[i + 1], xᵢ, c = :gray, alpha = 0.2)
    vline!(p[i + 1], xⱼ, c = :gray, alpha = 0.2)
end
plot(p..., size = (800, 600))
```

!!! note
    To check if a point set is a $(t,m,d)$-net, you can use the function `istmsnet` defined in the [tests file](https://github.com/SciML/QuasiMonteCarlo.jl/blob/2dce9905e564a85e1280115cc8af071674fc7d80/test/runtests.jl#L34) of this package.
    It uses the excellent [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl) package.
