# """
# ```julia
# struct SectionSample{T} <: SamplingAlgorithm
# ```
# `SectionSample(x0, sampler)` where `sampler` is any sampler above and `x0` is a vector of either `NaN` for a free dimension or some scalar for a constrained dimension.
# """
# struct SectionSample{T} <: SamplingAlgorithm
#     x0::Vector{T}
#     sa::SamplingAlgorithm
# end

# """
#     SectionSample{T}(x0, sa)

# `SectionSample(x0, sampler)` where `sampler` is any sampler above and `x0` is a vector of either `NaN` for a free dimension or some scalar for a constrained dimension.
# """
# struct SectionSample{R<:Real,I<:Integer,VR<:AbstractVector{R},VI<:AbstractVector{I}} <: SamplingAlgorithm
#     x0::VR
#     sa::SamplingAlgorithm
#     fixed_dims::VI
# end

# fixed_dimensions(section_sampler::SectionSample)::Vector{Int64} = findall(x -> x == false,
#                                                                           isnan.(section_sampler.x0))

# free_dimensions(section_sampler::SectionSample)::Vector{Int64} = findall(x -> x == true,
#                                                                          isnan.(section_sampler.x0))

# """
#     sample(n,lb,ub,K::SectionSample)

# Returns Tuples constrained to a section.
# In surrogate-based identification and control, optimization can alternate between unconstrained sampling in the full-dimensional parameter space, and sampling constrained on specific sections (e.g. a planes in a 3D volume),
# A SectionSample allows sampling and optimizing on a subset of 'free' dimensions while keeping 'fixed' ones constrained.
# The sampler is defined as in e.g.
# `section_sampler_y_is_10 = SectionSample([NaN64, NaN64, 10.0, 10.0], UniformSample())`
# where the first argument is a Vector{T} in which numbers are fixed coordinates and `NaN`s correspond to free dimensions, and the second argument is a SamplingAlgorithm which is used to sample in the free dimensions.
# """
# function sample(n::Integer, lb::Union{Number, Tuple, AbstractVector},
#                 ub::Union{Number, Tuple, AbstractVector}, section_sampler::SectionSample)
#     if n <= 0
#         throw(ZeroSamplesError())
#     end
#     if !check_bounds(lb, ub)
#         throw(UbLbWrong())
#     end
#     if lb isa Number
#         if isnan(section_sampler.x0[1])
#             return sample(n, lb, ub, section_sampler.sa)
#         else
#             return fill(section_sampler.x0[1], n)
#         end
#     else
#         d_free = free_dimensions(section_sampler)
#         new_samples = sample(n, lb[d_free], ub[d_free], section_sampler.sa)
#         out_as_vec = collect(repeat(section_sampler.x0', n, 1)')
#         for y in 1:size(out_as_vec, 2)
#             for (xi, d) in enumerate(d_free)
#                 out_as_vec[d, y] = new_samples[xi, y]
#             end
#         end
#         return out_as_vec
#     end
# end

# SectionSample(x0::AbstractVector, sa::SamplingAlgorithm) =
#     SectionSample(x0, sa, findall(isnan, x0))

# """
#     SectionSample(n, d, K::SectionSample)

# In surrogate-based identification and control, optimization can alternate between unconstrained sampling in the full-dimensional parameter space, and sampling constrained on specific sections (e.g. planes in a 3D volume).

# `SectionSample` allows sampling and optimizing on a subset of 'free' dimensions while keeping 'fixed' ones constrained.

# The sampler is defined

# `SectionSample([NaN64, NaN64, 10.0, 10.0], UniformSample())`

# where the first argument is a Vector{T} in which numbers are fixed coordinates and `NaN`s correspond to free dimensions, and the second argument is a SamplingAlgorithm which is used to sample in the free dimensions.
# """
# function sample(n::Integer, d::Integer, section_sampler::SectionSample, T=eltype(section_sampler.x0))
#     _check_sequence(n)
#     @assert eltype(section_sampler.x0) == T
#     @assert length(section_sampler.fixed_dims) == d
#     return sample(n, section_sampler)
# end

# @views function sample(n::Integer, section_sampler::SectionSample{T}) where T
#     samples = Matrix{T}(undef, n, length(section_sampler.x0))
#     fixed_dims = section_sampler.fixed_dims
#     samples[:,fixed_dims] .= sample(n, length(fixed_dims), section_sampler.sa, T)
#     return samples
# end
