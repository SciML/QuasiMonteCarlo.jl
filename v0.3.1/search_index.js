var documenterSearchIndex = {"docs":
[{"location":"design_matrix/#DesignMatrices","page":"Design Matrices","title":"Design Matrices","text":"","category":"section"},{"location":"design_matrix/#API","page":"Design Matrices","title":"API","text":"","category":"section"},{"location":"design_matrix/","page":"Design Matrices","title":"Design Matrices","text":"It is often convenient to generate multiple independent sequence, for error estimation (uncertainty quantification). The resulting sequences can be stored in what is often called a design matrix. In this package, this is achieved with the generate_design_matrices(n, d, ::DeterministicSamplingAlgorithm), ::RandomizationMethod, num_mats) function. num_mats is the number of independent realization. The resulting design matrix is a vector of matrix of length num_mats.","category":"page"},{"location":"design_matrix/","page":"Design Matrices","title":"Design Matrices","text":"QuasiMonteCarlo.generate_design_matrices","category":"page"},{"location":"design_matrix/#QuasiMonteCarlo.generate_design_matrices","page":"Design Matrices","title":"QuasiMonteCarlo.generate_design_matrices","text":"generate_design_matrices(n, d, sample_method::DeterministicSamplingAlgorithm,\nnum_mats, T = Float64)\ngenerate_design_matrices(n, d, sample_method::RandomSamplingAlgorithm,\nnum_mats, T = Float64)\ngenerate_design_matrices(n, lb, ub, sample_method,\nnum_mats = 2)\n\nCreate num_mats matrices each containing a QMC point set, where:\n\nn is the number of points to sample.\nd is the dimensionality of the point set in [0, 1)ᵈ,\nsample_method is the quasi-Monte Carlo sampling strategy used to create a deterministic point set out.\nT is the eltype of the point sets. For some QMC methods (Faure, Sobol) this can be Rational\n\nIf the bound lb and ub are specified instead of d, the samples will be transformed into the box [lb, ub].\n\n\n\n\n\ngenerate_design_matrices(n, d, sampler, R::NoRand, num_mats, T = Float64)\n\nR = NoRand() produces num_mats matrices each containing a different deterministic point set in [0, 1)ᵈ. Note that this is an ad hoc way to produce i.i.d sequence as it creates a deterministic point in dimension d × num_mats and split it in num_mats point set of dimension d.  This does not have any QMC garantuees.\n\n\n\n\n\n","category":"function"},{"location":"design_matrix/","page":"Design Matrices","title":"Design Matrices","text":"warning: Warning\nThe method generate_design_matrices(n, d, sampler, R::NoRand, num_mats, T = Float64) is an ad hoc way to produce a Design Matrix. Indeed, it creates a deterministic point set in dimension d × num_mats and splits it into num_mats point set of dimension d. The resulting sequences have no QMC guarantees. This seems to have been proposed in Section 5.1 of Saltelli, A. et al. (2010)[1] to do uncertainty quantification.","category":"page"},{"location":"design_matrix/","page":"Design Matrices","title":"Design Matrices","text":"[1]: Saltelli, A., Annoni, P., Azzini, I., Campolongo, F., Ratto, M., & Tarantola, S. (2010). Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index. Computer physics communications, 181(2), 259-270.","category":"page"},{"location":"design_matrix/#Example","page":"Design Matrices","title":"Example","text":"","category":"section"},{"location":"design_matrix/","page":"Design Matrices","title":"Design Matrices","text":"using QuasiMonteCarlo, Random\nRandom.seed!(1234)\nm = 4\nd = 3\nb = QuasiMonteCarlo.nextprime(d)\nN = b^m # Number of points\npad = m # Can also choose something as `2m` to get \"better\" randomization\n\n# 5 independent Randomized Faure sequences\nQuasiMonteCarlo.generate_design_matrices(N, d, FaureSample(R = OwenScramble(base = b, pad = m)), 5)","category":"page"},{"location":"randomization/#Randomization","page":"Randomization methods","title":"Randomization methods","text":"","category":"section"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"Most of the methods presented in Sampler are deterministic, i.e. X = sample(n, d, ::DeterministicSamplingAlgorithm) will always produce the same sequence X = (X_1 dots X_n).","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"The main issue with deterministic Quasi Monte Carlo sampling is that it does not allow easy error estimation as opposed to plain Monte Carlo where the variance can be estimated.","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"A Randomized Quasi Monte Carlo method must respect the two following criteria:","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"Have X_isim mathbbU(01^d) for each iin 1cdots n.\nPreserve the QMC properties, i.e. the randomized X still has low discrepancy.","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"This randomized version is unbiased and can be used to obtain confidence interval or to do sensitivity analysis.","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"A good reference is the book by A. Owen, especially the Chapters 15, 16 and 17.","category":"page"},{"location":"randomization/#API-for-randomization","page":"Randomization methods","title":"API for randomization","text":"","category":"section"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"abstract type RandomizationMethod end","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"There are two ways to obtain a randomized sequence:","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"Either directly use QuasiMonteCarlo.sample(n, d, DeterministicSamplingAlgorithm(R = SomeRandomizationMethod())) or sample(n, lb, up, DeterministicSamplingAlgorithm(R = RandomizationMethod())).\nOr, given n points d-dimensional points, all in 01^d one can do randomize(X, SomeRandomizationMethod()) where X is a dtimes n-matrix.","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"randomize","category":"page"},{"location":"randomization/#QuasiMonteCarlo.randomize","page":"Randomization methods","title":"QuasiMonteCarlo.randomize","text":"randomize(x, R::Shift)\n\nCranley Patterson Rotation i.e. y = (x .+ U) mod 1 where U ∼ 𝕌([0,1]ᵈ) and x is a d×n matrix\n\n\n\n\n\nrandomize(x, R::ScrambleMethod)\n\nReturn a scrambled version of x.  The scramble methods implemented are \n\nDigitalShift.\nOwenScramble: Nested Uniform Scramble which was introduced in Owen (1995).\nMatousekScramble: Linear Matrix Scramble which was introduced in Matousek (1998).\n\n\n\n\n\n","category":"function"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"The default method of DeterministicSamplingAlgorithm is NoRand","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"NoRand","category":"page"},{"location":"randomization/#QuasiMonteCarlo.NoRand","page":"Randomization methods","title":"QuasiMonteCarlo.NoRand","text":"NoRand <: RandomizationMethod\n\nNo Randomization is performed on the sampled sequence.\n\n\n\n\n\n","category":"type"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"To obtain multiple independent randomization of a sequence, i.e. Design Matrices, look at the Design Matrices section.","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"note: Note\nIn most other QMC packages, randomization is performed \"online\" as the points are samples. Here, randomization is performed after the deterministic sequence is generated. Both methods are useful in different contexts, the former is generally faster to produce one randomized sequence, while the latter is faster to produce independent realization of the sequence.PRs are welcomed to add \"online\" version of the sequence! See this comment for inspiration.Another way to view the two approaches is: given a computational budget of N points, one canPut all of it into a sequence of size, N, thus having the best estimator hatmu_N. The price to pay is that this estimation is not associated with a variance estimation.\nDivided your computational budget into N = ntimes M to get M independent estimator hatmu_n. From there one can compute an empirical variance of the estimator.","category":"page"},{"location":"randomization/#Scrambling-methods","page":"Randomization methods","title":"Scrambling methods","text":"","category":"section"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"abstract type ScrambleMethod <: RandomizationMethod end","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"ScrambleMethod","category":"page"},{"location":"randomization/#QuasiMonteCarlo.ScrambleMethod","page":"Randomization methods","title":"QuasiMonteCarlo.ScrambleMethod","text":"ScrambleMethod <: RandomizationMethod\n\nA scramble method needs at lease the scrambling base b, the number of \"bits\" to use pad (pad=32 is the default) and a seed rng (rng = Random.GLOBAL_RNG is the default). The scramble methods implementer are \n\nDigitalShift.\nOwenScramble: Nested Uniform Scramble which was introduced in Owen (1995).\nMatousekScramble: Linear Matrix Scramble which was introduced in Matousek (1998).\n\n\n\n\n\n","category":"type"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"ScramblingMethods(b, pad, rng) are well suited for (tmd)-nets in base b. b is the base used to scramble and pad the number of bits in b-ary decomposition, i.e. y simeq sum_k=1^textttpad y_ktextttb^k.","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"The pad is generally chosen as gtrsim log_b(n).","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"warning: Warning\nIn principle, the base b used for scrambling methods ScramblingMethods(b, pad, rng) can be an arbitrary integer. However, to preserve good Quasi Monte Carlo properties, it must match the base of the sequence to scramble. For example, (deterministic) Sobol' sequence are base b=2, (tmd) sequences while (deterministic) Faure sequences are (tmd) sequences in prime base i.e. b is an arbitrary prime number.","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"The implemented ScramblingMethods are","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"DigitalShift the simplest and faster method. For a point xin 01^d it does y_k = (x_k + U_k) mod b where U_k sim mathbbU(0 cdots b-1)","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"DigitalShift","category":"page"},{"location":"randomization/#QuasiMonteCarlo.DigitalShift","page":"Randomization methods","title":"QuasiMonteCarlo.DigitalShift","text":"DigitalShift <: ScrambleMethod\n\nDigital shift.  randomize(x, R::DigitalShift) returns a scrambled version of x. \n\nThe scramble method is Digital Shift. It scramble each corrdinate in base b as yₖ = (xₖ + Uₖ) mod b where Uₖ ∼ 𝕌({0:b-1}).  U is the same for every point points but i.i.d along every dimensions.\n\n\n\n\n\n","category":"type"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"MatousekScramble a.k.a. Linear Matrix Scramble is what people use in practice. Indeed, the observed performances are similar to OwenScramble for a lesser numerical cost.","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"MatousekScramble","category":"page"},{"location":"randomization/#QuasiMonteCarlo.MatousekScramble","page":"Randomization methods","title":"QuasiMonteCarlo.MatousekScramble","text":"MatousekScramble <: ScrambleMethod\n\nLinear Matrix Scramble aka Matousek' scramble.\n\nrandomize(x, R::MatousekScramble) returns a scrambled version of x.  The scramble method is Linear Matrix Scramble which was introduced in Matousek (1998). pad is the number of bits used for each points. One need pad ≥ log(base, n). \n\nReferences: Matoušek, J. (1998). On thel2-discrepancy for anchored boxes. Journal of Complexity, 14(4), 527-556.\n\n\n\n\n\n","category":"type"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"OwenScramble a.k.a. Nested Uniform Scramble is the most understood theoretically but is more costly to operate.","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"OwenScramble","category":"page"},{"location":"randomization/#QuasiMonteCarlo.OwenScramble","page":"Randomization methods","title":"QuasiMonteCarlo.OwenScramble","text":"OwenScramble <: ScrambleMethod\n\nNested Uniform Scramble aka Owen' scramble.\n\nrandomize(x, R::OwenScramble) returns a scrambled version of x.  The scramble method is Nested Uniform Scramble which was introduced in Owen (1995). pad is the number of bits used for each points. One needs pad ≥ log(base, n). \n\nReferences: Owen, A. B. (1995). Randomly permuted (t, m, s)-nets and (t, s)-sequences. In Monte Carlo and Quasi-Monte Carlo Methods in Scientific Computing: Proceedings of a conference at the University of Nevada, Las Vegas, Nevada, USA, June 23–25, 1994 (pp. 299-317). Springer New York.\n\n\n\n\n\n","category":"type"},{"location":"randomization/#Other-methods","page":"Randomization methods","title":"Other methods","text":"","category":"section"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"Shift(rng) a.k.a. Cranley-Patterson Rotation. It is by far the fastest method, it is used in LatticeRuleScramble for example.","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"Shift","category":"page"},{"location":"randomization/#QuasiMonteCarlo.Shift","page":"Randomization methods","title":"QuasiMonteCarlo.Shift","text":"Shifting(rng::AbstractRNG = Random.GLOBAL_RNG) <: RandomizationMethod\n\nCranley-Patterson rotation aka Shifting\n\nReferences: Cranley, R., & Patterson, T. N. (1976). Randomization of number theoretic methods for multiple integration. SIAM Journal on Numerical Analysis, 13(6), 904-914.\n\n\n\n\n\n","category":"type"},{"location":"randomization/#Example","page":"Randomization methods","title":"Example","text":"","category":"section"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"Randomization of a Faure sequence with various methods.","category":"page"},{"location":"randomization/#Generation","page":"Randomization methods","title":"Generation","text":"","category":"section"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"using QuasiMonteCarlo, Random\nRandom.seed!(1234)\nm = 4\nd = 3\nb = QuasiMonteCarlo.nextprime(d)\nN = b^m # Number of points\npad = m # Can also choose something as `2m` to get \"better\" randomization\n\n# Unrandomized low discrepancy sequence\nx_faure = QuasiMonteCarlo.sample(N, d, FaureSample())\n\n# Randomized version\nx_uniform = rand(d, N) # plain i.i.d. uniform\nx_shift = randomize(x_faure, Shift())\nx_nus = randomize(x_faure, OwenScramble(base = b, pad = pad)) # equivalent to sample(N, d, FaureSample(R = OwenScramble(base = b, pad = pad)))\nx_lms = randomize(x_faure, MatousekScramble(base = b, pad = pad))\nx_digital_shift = randomize(x_faure, DigitalShift(base = b, pad = pad))","category":"page"},{"location":"randomization/#Visualization-of-different-methods","page":"Randomization methods","title":"Visualization of different methods","text":"","category":"section"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"Plot the resulting sequences along dimensions 1 and 3.","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"using Plots\n# Setting I like for plotting\ndefault(fontfamily=\"Computer Modern\", linewidth=1, label=nothing, grid=true, framestyle=:default)\n\nd1 = 1\nd2 = 3 # you can try every combination of dimension (d1, d2)\nsequences = [x_uniform, x_faure, x_shift, x_digital_shift, x_lms, x_nus]\nnames = [\"Uniform\", \"Faure (deterministic)\", \"Shift\", \"Digital Shift\", \"Matousek Scramble\", \"Owen Scramble\"]\np = [plot(thickness_scaling=1.5, aspect_ratio=:equal) for i in sequences]\nfor (i, x) in enumerate(sequences)\n    scatter!(p[i], x[d1, :], x[d2, :], ms=1.5, c=1, grid=false)\n    title!(names[i])\n    xlims!(p[i], (0, 1))\n    ylims!(p[i], (0, 1))\n    yticks!(p[i], [0, 1])\n    xticks!(p[i], [0, 1])\n    hline!(p[i], range(0, 1, step=1 / 4), c=:gray, alpha=0.2)\n    vline!(p[i], range(0, 1, step=1 / 4), c=:gray, alpha=0.2)\n    hline!(p[i], range(0, 1, step=1 / 2), c=:gray, alpha=0.8)\n    vline!(p[i], range(0, 1, step=1 / 2), c=:gray, alpha=0.8)\nend\nplot(p..., size=(800, 600))","category":"page"},{"location":"randomization/#(t,m,d)-net-visualization","page":"Randomization methods","title":"(tmd)-net visualization","text":"","category":"section"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"Faure nets and its scrambled versions are digital (tmd)-net, it means that they have strong equipartition properties. On the following plot, we can (visually) verify that with Nested Uniform Scrambling, it also works with Linear Matrix Scrambling and Digital Shift. You must see one point per rectangle of volume 1b^m. Points on the \"left\" border of rectangles are included while those on the \"right\" are excluded. See Chapter 15.7 and Figure 15.10 for more details.","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"d1 = 1 \nd2 = 3 # you can try every combination of dimension (d1, d2)\nx = x_nus # Owen Scramble, you can try x_lms and x_digital_shift\np = [plot(thickness_scaling=1.5, aspect_ratio=:equal) for i in 0:m]\nfor i in 0:m\n    j = m - i\n    xᵢ = range(0, 1, step=1 / b^(i))\n    xⱼ = range(0, 1, step=1 / b^(j))\n    scatter!(p[i+1], x[d1, :], x[d2, :], ms=1.5, c=1, grid=false)\n    xlims!(p[i+1], (0, 1.01))\n    ylims!(p[i+1], (0, 1.01))\n    yticks!(p[i+1], [0, 1])\n    xticks!(p[i+1], [0, 1])\n    hline!(p[i+1], xᵢ, c=:gray, alpha=0.2)\n    vline!(p[i+1], xⱼ, c=:gray, alpha=0.2)\nend\nplot(p..., size=(800, 600))","category":"page"},{"location":"randomization/","page":"Randomization methods","title":"Randomization methods","text":"note: Note\nTo check if a point set is a (tmd)-net, you can use the function istmsnet defined in the tests file of this package. It uses the excellent IntervalArithmetic.jl package.","category":"page"},{"location":"#QuasiMonteCarlo.jl:-Quasi-Monte-Carlo-(QMC)-Samples-Made-Easy","page":"Home","title":"QuasiMonteCarlo.jl: Quasi-Monte Carlo (QMC) Samples Made Easy","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"QuasiMonteCarlo.jl is a lightweight package for generating Quasi-Monte Carlo (QMC) samples using various different methods.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install QuasiMonteCarlo.jl, use the Julia package manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"QuasiMonteCarlo\")","category":"page"},{"location":"#Get-Started","page":"Home","title":"Get Started","text":"","category":"section"},{"location":"#Basic-API","page":"Home","title":"Basic API","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using QuasiMonteCarlo, Distributions\nlb = [0.1,-0.5]\nub = [1.0,20.0]\nn = 5\nd = 2\n\ns = QuasiMonteCarlo.sample(n,lb,ub,GridSample())\ns = QuasiMonteCarlo.sample(n,lb,ub,Uniform())\ns = QuasiMonteCarlo.sample(n,lb,ub,SobolSample())\ns = QuasiMonteCarlo.sample(n,lb,ub,LatinHypercubeSample())\ns = QuasiMonteCarlo.sample(n,lb,ub,LatticeRuleSample())\ns = QuasiMonteCarlo.sample(n,lb,ub,HaltonSample())","category":"page"},{"location":"","page":"Home","title":"Home","text":"The output s is a matrix, so one can use things like @uview from UnsafeArrays.jl for a stack-allocated view of the ith point:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using UnsafeArrays\n@uview s[:,i]","category":"page"},{"location":"#MC-vs-QMC","page":"Home","title":"MC vs QMC","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"We illustrate the gain of QMC methods over plain Monte Carlo using the 5-dimensional example from Section 15.9 in the book by A. Owen.","category":"page"},{"location":"","page":"Home","title":"Home","text":"f₁(𝐱) = prod(1 + √(12)/5*(xⱼ - 1/2) for xⱼ ∈ 𝐱)\nμ_exact = 1 # = ∫ f₁(𝐱) d⁵𝐱","category":"page"},{"location":"","page":"Home","title":"Home","text":"One can estimate the integral mu using plain Monte Carlo, or Quasi Monte Carlo or Randomized Quasi Monte Carlo. See the other section of this documentation for more information on the functions used in the example.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using QuasiMonteCarlo, Random, Distributions\nusing Plots; default(fontfamily=\"Computer Modern\")\nRandom.seed!(1234)\nd = 5 # Dimension (= prime base for Faure net)\nb = 2 # Base for Sobol net\nm_max = 19\nm_max_Faure = 8\nN = b^m_max\n\n# Generate sequences\nseq_MC = QuasiMonteCarlo.sample(N, d, Uniform()) # Monte Carlo i.i.d Uniform sampling\nseq_QMC_Sobol = QuasiMonteCarlo.sample(N, d, SobolSample()) # Sobol net\nseq_RQMC_Sobol = QuasiMonteCarlo.sample(N, d, SobolSample(R = OwenScramble(base = b, pad = 32))) # Randomized version of Sobol net\nseq_RQMC_Faure = QuasiMonteCarlo.sample(d^m_max_Faure, d, FaureSample(R = OwenScramble(base = d, pad = 32))) # Randomized version of Faure net\n\n# Estimate the integral for different n with different estimator μ̂ₙ\nμ_MC = [mean(f₁(x) for x in eachcol(seq_MC[:, 1:b^m])) for m in 1:m_max]\nμ_QMC_Sobol = [mean(f₁(x) for x in eachcol(seq_QMC_Sobol[:, 1:b^m])) for m in 1:m_max]\nμ_RQMC_Sobol = [mean(f₁(x) for x in eachcol(seq_RQMC_Sobol[:, 1:b^m])) for m in 1:m_max]\nμ_RQMC_Faure = [mean(f₁(x) for x in eachcol(seq_RQMC_Faure[:, 1:d^m])) for m in 1:m_max_Faure]\n\n# Plot the error |μ̂-μ| vs n\nplot(b.^(1:m_max), abs.(μ_MC .- μ_exact), label=\"MC\")\nplot!(b.^(1:m_max), abs.(μ_QMC_Sobol .- μ_exact), label=\"QMC Sobol\")\nplot!(b.^(1:m_max), abs.(μ_RQMC_Sobol .- μ_exact), label=\"RQMC Sobol\")\nplot!(d .^(1:m_max_Faure), abs.(μ_RQMC_Faure .- μ_exact), label=\"RQMC Faure\")\nplot!(n -> n^(-1/2), b.^(1:m_max), c = :black, s = :dot, label = \"n^(-1/2)\")\nplot!(n -> n^(-3/2), b.^(1:m_max), c = :black, s = :dash, label = \"n^(-3/2)\") \n# n^(-3/2) is the theoretical scaling for scrambled nets e.g. Theorem 17.5. in https://artowen.su.domains/mc/qmcstuff.pdf\nxlims!(1, 1e6)\nylims!(1e-9, 1)\nxaxis!(:log10)\nyaxis!(:log10)\nxlabel!(\"n\", legend = :bottomleft)\nylabel!(\"|μ̂-μ|\")","category":"page"},{"location":"#Adding-a-new-sampling-method","page":"Home","title":"Adding a new sampling method","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Adding a new sampling method is a two-step process:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Add a new SamplingAlgorithm type.\nOverload the sample function with the new type.","category":"page"},{"location":"","page":"Home","title":"Home","text":"All sampling methods are expected to return a matrix with dimension d by n, where d is the dimension of the sample space and n is the number of samples.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Example","category":"page"},{"location":"","page":"Home","title":"Home","text":"struct NewAmazingSamplingAlgorithm{OPTIONAL} <: SamplingAlgorithm end\n\nfunction sample(n,lb,ub,::NewAmazingSamplingAlgorithm)\n    if lb isa Number\n        ...\n        return x\n    else\n        ...\n        return reduce(hcat, x)\n    end\nend","category":"page"},{"location":"#Contributing","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Please refer to the SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages for guidance on PRs, issues, and other matters relating to contributing to SciML.\nThere are a few community forums:\nthe #diffeq-bridged channel in the Julia Slack\nJuliaDiffEq on Gitter\non the Julia Discourse forums\nsee also SciML Community page","category":"page"},{"location":"#Reproducibility","page":"Home","title":"Reproducibility","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg # hide\nPkg.status() # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"<details><summary>and using this machine and Julia version.</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using InteractiveUtils # hide\nversioninfo() # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg # hide\nPkg.status(;mode = PKGMODE_MANIFEST) # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can also download the\n<a href=\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"using TOML\nversion = TOML.parse(read(\"../../Project.toml\",String))[\"version\"]\nname = TOML.parse(read(\"../../Project.toml\",String))[\"name\"]\nlink = \"https://github.com/SciML/\"*name*\".jl/tree/gh-pages/v\"*version*\"/assets/Manifest.toml\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"\">manifest</a> file and the\n<a href=\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"using TOML\nversion = TOML.parse(read(\"../../Project.toml\",String))[\"version\"]\nname = TOML.parse(read(\"../../Project.toml\",String))[\"name\"]\nlink = \"https://github.com/SciML/\"*name*\".jl/tree/gh-pages/v\"*version*\"/assets/Project.toml\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"\">project</a> file.","category":"page"},{"location":"samplers/#Sampler-APIs","page":"Sampler APIs","title":"Sampler APIs","text":"","category":"section"},{"location":"samplers/#Sample","page":"Sampler APIs","title":"Sample","text":"","category":"section"},{"location":"samplers/","page":"Sampler APIs","title":"Sampler APIs","text":"QuasiMonteCarlo.sample","category":"page"},{"location":"samplers/#QuasiMonteCarlo.sample","page":"Sampler APIs","title":"QuasiMonteCarlo.sample","text":"sample(n::Integer, d::Integer, S::SamplingAlgorithm, T = Float64)\nsample(n::Integer, lb::T, ub::T, S::SamplingAlgorithm) where T <: Union{Base.AbstractVecOrTuple, Number}\n\nReturn a QMC point set where:\n\nn is the number of points to sample.\nS is the quasi-Monte Carlo sampling strategy. \n\nThe point set is in a d-dimensional unit box [0, 1]^d.  If the bounds are specified, the sample is transformed (translation + scaling) into a box [lb, ub] where:\n\nlb is the lower bound for each variable. Its length fixes the dimensionality of the sample.\nub is the upper bound. Its dimension must match length(lb).\n\nIn the first method the type of the point set is specified by T while in the second method the output type is infered from the bound types.\n\n\n\n\n\nsample(n::Integer, lb::T, ub::T, D::Distributions.Sampleable, T = eltype(D))\nsample(n::Integer, lb::T, ub::T, D::Distributions.Sampleable) where T <: Union{Base.AbstractVecOrTuple, Number}\n\nReturn a point set from a distribution D:\n\nn is the number of points to sample.\nD is a Distributions.Sampleable from Distributions.jl.\n\nThe point set is in a d-dimensional unit box [0, 1]^d.  If the bounds are specified instead of just d, the sample is transformed (translation + scaling) into a box [lb, ub] where:\n\nlb is the lower bound for each variable. Its length fixes the dimensionality of the sample.\nub is the upper bound. Its dimension must match length(lb).\n\n\n\n\n\n","category":"function"},{"location":"samplers/#Samplers","page":"Sampler APIs","title":"Samplers","text":"","category":"section"},{"location":"samplers/","page":"Sampler APIs","title":"Sampler APIs","text":"Samplers are divided into two subtypes","category":"page"},{"location":"samplers/","page":"Sampler APIs","title":"Sampler APIs","text":"abstract type SamplingAlgorithm end\nabstract type RandomSamplingAlgorithm <: SamplingAlgorithm end\nabstract type DeterministicSamplingAlgorithm <: SamplingAlgorithm end","category":"page"},{"location":"samplers/#Deterministic-Sampling-Algorithm","page":"Sampler APIs","title":"Deterministic Sampling Algorithm","text":"","category":"section"},{"location":"samplers/","page":"Sampler APIs","title":"Sampler APIs","text":"All DeterministicSamplingAlgorithm have NoRand() as their default RandomizationMethod, see Randomization methods and Design Matrices section for more information on randomization.","category":"page"},{"location":"samplers/","page":"Sampler APIs","title":"Sampler APIs","text":"GridSample","category":"page"},{"location":"samplers/#QuasiMonteCarlo.GridSample","page":"Sampler APIs","title":"QuasiMonteCarlo.GridSample","text":"GridSample(R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm\n\nA simple rectangular grid lattice.\n\nIn more than 2 dimensions, grids have worse discrepancy than simple Monte Carlo. As a result, they should almost never be used for multivariate integration; their use is as a starting point for other algorithms.\n\n\n\n\n\n","category":"type"},{"location":"samplers/","page":"Sampler APIs","title":"Sampler APIs","text":"SobolSample","category":"page"},{"location":"samplers/#QuasiMonteCarlo.SobolSample","page":"Sampler APIs","title":"QuasiMonteCarlo.SobolSample","text":"SobolSample(R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm\n\nSamples taken from Sobol's base-2 sequence.\n\n\n\n\n\n","category":"type"},{"location":"samplers/","page":"Sampler APIs","title":"Sampler APIs","text":"warning: Warning\nThe QuasiMonteCarlo.jl package relies on the Sobol.jl package to sample Sobol nets. The choice, there is to NOT start the sequence at 0. This is debatable, see this issue and ref therein for more context.","category":"page"},{"location":"samplers/","page":"Sampler APIs","title":"Sampler APIs","text":"FaureSample\nLatticeRuleSample\nHaltonSample\nGoldenSample\nKroneckerSample","category":"page"},{"location":"samplers/#QuasiMonteCarlo.FaureSample","page":"Sampler APIs","title":"QuasiMonteCarlo.FaureSample","text":"FaureSample(R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm\n\nA Faure low-discrepancy sequence.\n\nFaure-distributed samples cover all dimensions evenly, using the same set of points for all variables, up to ordering.\n\nWhen scrambled, randomized Faure sequences provide worst-case guarantees that variance will be at most exp(1) ≈ 2.718 times greater than for a purely Monte Carlo integral. However, they are much less efficient than the Sobol sequence at integrating functions with low effective dimension (functions where the first few inputs dominate the evaluation).\n\nThe Faure sequence in dimension s forms a (0, s)-sequence with base b = nextprime(s).\n\nA Faure sequence must have length k * base^s with k < base < 1.\n\nReferences: Faure, H. (1982). Discrépance de suites associées à un système de numération (en dimension s). Acta Arith., 41, 337-351. Owen, A. B. (1997). Monte Carlo variance of scrambled net quadrature. SIAM Journal on Numerical Analysis, 34(5), 1884-1910.\n\n\n\n\n\n","category":"type"},{"location":"samplers/#QuasiMonteCarlo.LatticeRuleSample","page":"Sampler APIs","title":"QuasiMonteCarlo.LatticeRuleSample","text":"LatticeRuleSample() <: DeterministicSamplingAlgorithm\n\nGenerate a point set using a lattice rule.\n\n\n\n\n\n","category":"type"},{"location":"samplers/#QuasiMonteCarlo.HaltonSample","page":"Sampler APIs","title":"QuasiMonteCarlo.HaltonSample","text":"HaltonSample(R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm\n\nCreate a Halton sequence.\n\n\n\n\n\n","category":"type"},{"location":"samplers/#QuasiMonteCarlo.GoldenSample","page":"Sampler APIs","title":"QuasiMonteCarlo.GoldenSample","text":"GoldenSample()\n\nGenerate a quasirandom Kronecker sequence using powers of the generalized golden ratio.\n\nThe harmonious, or generalized golden, ratios are defined as the solutions to the equation: x^d = x + 1\n\nWhere d is the dimension of the sequence. The Golden sequence is then equivalent to Kronecker([x^-i for i in 1:d]).\n\nWARNING: the generalized golden sequence in more than 2 dimensions is purely experimental. It likely has poor discrepancy in high dimensions, and should not be used without verifying answers against a better-known quasirandom sequence. Try a rank-1 lattice rule instead.\n\nReferences: Roberts, M. (2018). The Unreasonable Effectiveness of Quasirandom Sequences. Extreme Learning. http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/\n\n\n\n\n\n","category":"function"},{"location":"samplers/#QuasiMonteCarlo.KroneckerSample","page":"Sampler APIs","title":"QuasiMonteCarlo.KroneckerSample","text":"KroneckerSample(generator::AbstractVector, R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm\n\nA Kronecker sequence is a point set generated using a vector and equation: x[i] = i * generator .% 1\n\nWhere i runs from 1 through the sample size n. This sequence will be equidistributed (uniform in the infinite limit) so long as the components of generator are linearly independent over the field of rational numbers.\n\nIf no generator is specified, a lattice based on the generalized golden ratio is used; see GoldenSample for more information.\n\nKronecker sequences are not recommended for use in more than 3 dimensions, as theory on them is sparse. LatticeRuleSample will return rank-1 lattice rules, which behave similarly to Kronecker sequences but have better properties.\n\nReferences: Leobacher, G., & Pillichshammer, F. (2014). Introduction to quasi-Monte Carlo integration and applications. Switzerland: Springer International Publishing. https://link.springer.com/content/pdf/10.1007/978-3-319-03425-6.pdf\n\n\n\n\n\n","category":"type"},{"location":"samplers/#Random-Sampling-Algorithm","page":"Sampler APIs","title":"Random Sampling Algorithm","text":"","category":"section"},{"location":"samplers/","page":"Sampler APIs","title":"Sampler APIs","text":"LatinHypercubeSample","category":"page"},{"location":"samplers/#QuasiMonteCarlo.LatinHypercubeSample","page":"Sampler APIs","title":"QuasiMonteCarlo.LatinHypercubeSample","text":"LatinHypercubeSample(rng::AbstractRNG = Random.GLOBAL_RNG) <: RandomSamplingAlgorithm\n\nA Latin Hypercube is a point set with the property that every one-dimensional interval (i/n, i+1/n) contains exactly one point. This is a good way to sample a high-dimensional space, as it is more uniform than a random sample but does not require as many points as a full net.\n\n\n\n\n\n","category":"type"}]
}
