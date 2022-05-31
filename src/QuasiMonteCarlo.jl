module QuasiMonteCarlo

using Sobol, LatinHypercubeSampling, LatticeRules, Distributions

abstract type SamplingAlgorithm end

"""
GridSample{T}
T is the step dx for lb:dx:ub
"""
struct GridSample{T} <: SamplingAlgorithm
    dx::T
end

struct UniformSample <: SamplingAlgorithm end
struct SobolSample <: SamplingAlgorithm end
struct LatinHypercubeSample <: SamplingAlgorithm end
struct LatticeRuleSample <: SamplingAlgorithm end

"""
LowDiscrepancySample{T}
T is the base for the sequence
"""
struct LowDiscrepancySample{T} <: SamplingAlgorithm
    base::T
end

struct RandomSample <: SamplingAlgorithm end

struct KroneckerSample{A,B} <: SamplingAlgorithm
    alpha::A
    s0::B
end

struct GoldenSample <: SamplingAlgorithm end

struct SectionSample{T} <: SamplingAlgorithm
    x0::Vector{T}
    sa::QuasiMonteCarlo.SamplingAlgorithm
end

"""
sample(n,lb,ub,S::GridSample)
Returns a tuple containing numbers in a grid.
"""
function sample(n,lb,ub,S::GridSample)
    dx = S.dx
    if lb isa Number
        return vec(rand(lb:S.dx:ub,n))
    else
        d = length(lb)
        x = [[rand(lb[j]:dx[j]:ub[j]) for j = 1:d] for i in 1:n]
        return reduce(hcat,x)
    end
end

"""
sample(n,lb,ub,::UniformRandom)
Returns a tuple containing uniform random numbers.
"""
function sample(n,lb,ub,::UniformSample)
    if lb isa Number
        return rand(Uniform(lb,ub),n)
    else
        d = length(lb)
        x = [[rand(Uniform(lb[j],ub[j])) for j in 1:d] for i in 1:n]
        return reduce(hcat,x)
    end
end

"""
sample(n,lb,ub,::SobolSampling)
Returns a tuple containing Sobol sequences.
"""
function sample(n,lb,ub,::SobolSample)
    s = SobolSeq(lb,ub)
    skip(s,n)
    if lb isa Number
        return [next!(s)[1] for i = 1:n]
    else
        return reduce(hcat,[next!(s) for i = 1:n])
    end
end

"""
sample(n,lb,ub,::LatinHypercube)
Returns a tuple containing LatinHypercube sequences.
"""
function sample(n,lb,ub,::LatinHypercubeSample)
    d = length(lb)
    if lb isa Number
        x = vec(LHCoptim(n,d,1)[1])
        # x∈[0,n], so affine transform
        return @. (ub-lb) * x/(n) + lb
    else
        lib_out = collect(float(LHCoptim(n,d,1)[1])')
        # x∈[0,n], so affine transform column-wise
        @inbounds for c = 1:d
            lib_out[c, :] = (ub[c]-lb[c])*lib_out[c, :]/n .+ lb[c]
        end
        return lib_out
    end
end

"""
sample(n,lb,ub,::LatticeRuleSample)
Returns a matrix with the `n` rank-1 lattice points in each column if `lb` is a vector, or a vector with the `n` rank-1 lattice points if `lb` is a number.
"""
function sample(n,lb,ub,::LatticeRuleSample)
    if lb isa Number
        lat = ShiftedLatticeRule(1)
        pts = reduce(vcat,Iterators.take(lat, n))
        # transform from [0, 1] to [lb, ub]
        @inbounds for i in 1:n
            pts[i] = (ub-lb)*pts[i] + lb
        end
        return pts
    else
        d = length(lb)
        lat = ShiftedLatticeRule(d)
        pts = reduce(hcat,Iterators.take(lat, n))
        # transform from [0, 1]^d to [lb, ub]
        @inbounds for j in 1:n
            for i in 1:d
                pts[i, j] = (ub[i]-lb[i])*pts[i, j] + lb[i]
            end
        end
        return pts
    end
end

"""
sample(n,lb,ub,S::LowDiscrepancySample)
Low-discrepancy sample:
- Dimension 1: Van der Corput sequence
- Dimension > 1: Halton sequence
If dimension d > 1, all bases must be coprime with one other.
"""
function sample(n,lb,ub,S::LowDiscrepancySample)
    @assert length(lb) == length(ub)

    d = length(lb)
    t = float(eltype(lb)) # infer data type for return values
    if d == 1
        #Van der Corput
        b = S.base
        x = zeros(t, n)
        for i = 1:n
            expansion = digits(i,base = b)
            L = length(expansion)
            val = zero(t)
            for k = 1:L
                val += expansion[k]*float(b)^(-(k-1)-1)
            end
            x[i] = val
        end
        # It is always defined on the unit interval, resizing:
        return @. (ub-lb) * x + lb
    else
        #Halton sequence
        x = zeros(t,d,n)
        for j = 1:d
            b = S.base[j]
            for i = 1:n
                val = zero(t)
                expansion = digits(i, base = b)
                L = length(expansion)
                for k = 1:L
                    val += expansion[k]*float(b)^(-(k-1)-1)
                end
                x[j,i] = val
            end
        end
        #Resizing
        # x∈[0,1], so affine transform column-wise
        @inbounds for c = 1:d
            x[c, :] = (ub[c]-lb[c])*x[c, :] .+ lb[c]
        end
        return x
    end
end

"""
sample(n,d,K::KroneckerSample)

Returns a Tuple containing numbers following the Kronecker sample
"""
function sample(n,lb,ub,K::KroneckerSample)
    d = length(lb)
    alpha = K.alpha
    s0 = K.s0
    if d == 1
        x = zeros(n)
        @inbounds for i = 1:n
            x[i] = (s0+i*alpha)%1
        end
        return @. (ub-lb) * x + lb
    else
        x = zeros(n,d)
        @inbounds for j = 1:d
            for i = 1:n
                x[i,j] = (s0[j] + i*alpha[j])%i
            end
        end
        #Resizing
        # x∈[0,1], so affine transform column-wise
        @inbounds for c = 1:d
            x[:,c] = (ub[c]-lb[c])*x[:,c] .+ lb[c]
        end

        y = collect(x')
        return y
    end
end

function sample(n,lb,ub, G::GoldenSample)
    d = length(lb)
    if d == 1
        x = zeros(n)
        g = (sqrt(5)+1)/2
        a = 1.0/g
        for i = 1:n
            x[i] = (0.5+a*i)%1
        end
        return @. (ub-lb) * x + lb
    else
        x = zeros(n,d)
        for j = 1:d
            #Approximate solution of x^(d+1) = x + 1, a simple newton is good enough
            y = 2.0
            for s = 1:10
                g = (1+y)^(1/(j+1))
            end
            a = 1.0/g
            for i = 1:n
                x[i,j] = (0.5+a*i)%1
            end
        end
        @inbounds for c = 1:d
            x[:,c] = (ub[c]-lb[c])*x[:,c] .+ lb[c]
        end
        y = collect(x')
        return y
    end
end

fixed_dimensions(
        section_sampler::SectionSample)::Vector{Int64} = findall(
    x->x == false, isnan.(section_sampler.x0))

free_dimensions(
        section_sampler::SectionSample)::Vector{Int64} = findall(
    x->x == true, isnan.(section_sampler.x0))

"""
sample(n,lb,ub,K::SectionSample)

Returns Tuples constrained to a section.

In surrogate-based identification and control, optimization can alternate between unconstrained sampling in the full-dimensional parameter space, and sampling constrained on specific sections (e.g. a planes in a 3D volume),

A SectionSampler allows sampling and optimizing on a subset of 'free' dimensions while keeping 'fixed' ones constrained.
The sampler is defined as in e.g.

`section_sampler_y_is_10 = SectionSample([NaN64, NaN64, 10.0, 10.0], UniformSample())`

where the first argument is a Vector{T} in which numbers are fixed coordinates and `NaN`s correspond to free dimensions, and the second argument is a SamplingAlgorithm which is used to sample in the free dimensions.
"""
function sample(n,lb,ub,section_sampler::SectionSample)
    if lb isa Number
        if isnan(section_sampler.x0[1])
            return sample(n, lb, ub, section_sampler.sa)
        else
            return fill(section_sampler.x0[1], n)
        end
    else
        d_free = free_dimensions(section_sampler)
        new_samples = sample(n, lb[d_free], ub[d_free], section_sampler.sa)
        out_as_vec = collect(repeat(section_sampler.x0', n, 1)')
        for y in 1:size(out_as_vec,2)
            for (xi, d) in enumerate(d_free)
                out_as_vec[d,y] = new_samples[y][xi]
            end
        end
        return out_as_vec
    end
end

"""
sample(n,d,D::Distribution)
Returns a tuple containing numbers distributed as D.
"""
function sample(n,d,D::Distribution)
    if d == 1
        return rand(D,n)
    else
        x = [[rand(D) for j in 1:d] for i in 1:n]
        return reduce(hcat, x)
    end
end

function generate_design_matrices(n,lb,ub,sampler,num_mats = 2)
    @assert length(lb) == length(ub)
    d = length(lb)
    out = sample(n, repeat(lb,num_mats), repeat(ub,num_mats),sampler)
    [out[(j*d+1):((j+1)*d),:] for j in 0:num_mats-1]
end

export GridSample, UniformSample, SobolSample, LatinHypercubeSample, LatticeRuleSample,
       RandomSample, LowDiscrepancySample, GoldenSample, KroneckerSample, SectionSample

end # module
