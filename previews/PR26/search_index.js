var documenterSearchIndex = {"docs":
[{"location":"#QuasiMonteCarlo.jl:-Quasi-Monte-Carlo-(QMC)-Samples-Made-Easy","page":"Home","title":"QuasiMonteCarlo.jl: Quasi-Monte Carlo (QMC) Samples Made Easy","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"QuasiMonteCarlo.jl is a lightweight package for generating Quasi-Monte Carlo (QMC) samples using various different methods.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install QuasiMonteCarlo.jl, use the Julia package manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"QuasiMonteCarlo\")","category":"page"},{"location":"#Example","page":"Home","title":"Example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using QuasiMonteCarlo, Distributions\nlb = [0.1,-0.5]\nub = [1.0,20.0]\nn = 5\nd = 2\n\ns = QuasiMonteCarlo.sample(n,lb,ub,GridSample([0.1,0.5]))\ns = QuasiMonteCarlo.sample(n,lb,ub,UniformSample())\ns = QuasiMonteCarlo.sample(n,lb,ub,SobolSample())\ns = QuasiMonteCarlo.sample(n,lb,ub,LatinHypercubeSample())\ns = QuasiMonteCarlo.sample(n,lb,ub,LatticeRuleSample())\ns = QuasiMonteCarlo.sample(n,lb,ub,LowDiscrepancySample([10,3]))","category":"page"},{"location":"","page":"Home","title":"Home","text":"The output s is a matrix, so one can use things like @uview from UnsafeArrays.jl for a stack-allocated view of the ith point:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using UnsafeArrays\n@uview s[:,i]","category":"page"},{"location":"#Adding-a-new-sampling-method","page":"Home","title":"Adding a new sampling method","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Adding a new sampling method is a two-step process:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Add a new SamplingAlgorithm type.\nOverload the sample function with the new type.","category":"page"},{"location":"","page":"Home","title":"Home","text":"All sampling methods are expected to return a matrix with dimension d by n, where d is the dimension of the sample space and n is the number of samples.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Example","category":"page"},{"location":"","page":"Home","title":"Home","text":"struct NewAmazingSamplingAlgorithm{OPTIONAL} <: SamplingAlgorithm end\n\nfunction sample(n,lb,ub,::NewAmazingSamplingAlgorithm)\n    if lb isa Number\n        ...\n        return x\n    else\n        ...\n        return reduce(hcat, x)\n    end\nend","category":"page"},{"location":"#Contributing","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Please refer to the SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages for guidance on PRs, issues, and other matters relating to contributing to SciML.\nThere are a few community forums:\nthe #diffeq-bridged channel in the Julia Slack\nJuliaDiffEq on Gitter\non the Julia Discourse forums\nsee also SciML Community page","category":"page"},{"location":"samplers/#Sampler-APIs","page":"Sampler APIs","title":"Sampler APIs","text":"","category":"section"},{"location":"samplers/#Sample","page":"Sampler APIs","title":"Sample","text":"","category":"section"},{"location":"samplers/","page":"Sampler APIs","title":"Sampler APIs","text":"sample","category":"page"},{"location":"samplers/#Samplers","page":"Sampler APIs","title":"Samplers","text":"","category":"section"},{"location":"samplers/","page":"Sampler APIs","title":"Sampler APIs","text":"GridSample\nUniformSample\nSobolSample\nLatinHypercubeSample\nLatticeRuleSample\nLowDiscrepancySample\nGoldenSample\nKroneckerSample\nSectionSample","category":"page"},{"location":"samplers/#QuasiMonteCarlo.GridSample","page":"Sampler APIs","title":"QuasiMonteCarlo.GridSample","text":"GridSample{T}\n\nThe grid is given by lb:dx[i]:ub in the ith direction.\n\n\n\n\n\n","category":"type"},{"location":"samplers/#QuasiMonteCarlo.UniformSample","page":"Sampler APIs","title":"QuasiMonteCarlo.UniformSample","text":"struct UniformSample <: SamplingAlgorithm end\n\nUniformly distributed random numbers.\n\n\n\n\n\n","category":"type"},{"location":"samplers/#QuasiMonteCarlo.SobolSample","page":"Sampler APIs","title":"QuasiMonteCarlo.SobolSample","text":"struct SobolSample <: SamplingAlgorithm end\n\nSamples using the Sobol sequence using Sobol.jl.\n\n\n\n\n\n","category":"type"},{"location":"samplers/#QuasiMonteCarlo.LatinHypercubeSample","page":"Sampler APIs","title":"QuasiMonteCarlo.LatinHypercubeSample","text":"struct LatinHypercubeSample <: SamplingAlgorithm end\n\nSamples using a Latin Hypercube using LatinHypercubeSampling.jl\n\nLatinHypercubeSample(threading=false)\n\nKeyword arguments:\n\nthreading: whether to use threading. Default is false, i.e. serial.\n\n\n\n\n\n","category":"type"},{"location":"samplers/#QuasiMonteCarlo.LatticeRuleSample","page":"Sampler APIs","title":"QuasiMonteCarlo.LatticeRuleSample","text":"struct LatticeRuleSample <: SamplingAlgorithm end\n\nSamples using a randomly-shifted rank-1 lattice rule using LatticeRules.jl\n\n\n\n\n\n","category":"type"},{"location":"samplers/#QuasiMonteCarlo.LowDiscrepancySample","page":"Sampler APIs","title":"QuasiMonteCarlo.LowDiscrepancySample","text":"struct LowDiscrepancySample{T} <: SamplingAlgorithm\n\nbase[i] is the base in the ith direction.\n\n\n\n\n\n","category":"type"},{"location":"samplers/#QuasiMonteCarlo.GoldenSample","page":"Sampler APIs","title":"QuasiMonteCarlo.GoldenSample","text":"struct GoldenSample <: SamplingAlgorithm end\n\n\n\n\n\n","category":"type"},{"location":"samplers/#QuasiMonteCarlo.KroneckerSample","page":"Sampler APIs","title":"QuasiMonteCarlo.KroneckerSample","text":"struct KroneckerSample{A,B} <: SamplingAlgorithm\n\nKroneckerSample(alpha, s0) for a Kronecker sequence, where alpha is an length-d vector of irrational numbers (often sqrt(d)) and s0 is a length-d seed vector (often 0).\n\n\n\n\n\n","category":"type"},{"location":"samplers/#QuasiMonteCarlo.SectionSample","page":"Sampler APIs","title":"QuasiMonteCarlo.SectionSample","text":"struct SectionSample{T} <: SamplingAlgorithm\n\nSectionSample(x0, sampler) where sampler is any sampler above and x0 is a vector of either NaN for a free dimension or some scalar for a constrained dimension.\n\n\n\n\n\n","category":"type"}]
}
