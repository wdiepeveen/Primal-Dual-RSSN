#
#      Rn - The manifold of the n-dimensional (real valued) Euclidean space
#
# Manopt.jl, R. Bergmann, 2019
import LinearAlgebra: I, norm
import Base: exp, log, show
export Euclidean, RnPoint, RnTVector
export distance, exp, log, norm, dot, manifoldDimension, show, getValue
export zeroTVector, tangentONB, randomMPoint, randomTVector, injectivityRadius
export project
export validateMPoint, validateTVector, typeofMPoint, typeofTVector
# Types
# ---

@doc doc"""
    Euclidean <: Manifold

The manifold $\mathcal M = \mathbb R^n$ of the $n$-dimensional Euclidean vector
space. We employ the notation $\langle\cdot,\cdot\rangle$ for the inner product
and $\lVert\cdot\rVert_2$ for its induced norm.

# Abbreviation

`Rn`

# Constructor

    Euclidean(n)

construct the n-dimensional Euclidean space $\mathbb R^n$.
"""
struct Euclidean <: Manifold
  name::String
  dimension::Int
  abbreviation::String
  Euclidean(dimension::Int) = new("$dimension-dimensional Euclidean space",dimension,"R$dimension")
end

@doc doc"""
    RnPoint <: MPoint

the point $x\in\mathcal M$ for $\mathcal M=\mathbb R^n$ represented by an
$n$-dimensional `Vector{T}`, where `T <: AbstractFloat`.
"""
struct RnPoint{T<:AbstractFloat} <: MPoint
  value::Vector{T}
  RnPoint{T}( value::Vector{T} ) where T<:AbstractFloat = new(value)
  RnPoint{T}( value::T ) where T<:AbstractFloat = new([value])
end
RnPoint(value::T) where {T <: AbstractFloat} = RnPoint{T}(value)
RnPoint(value::Vector{T}) where {T <: AbstractFloat} = RnPoint{T}(value)
getValue(x::RnPoint) = length(x.value)==1 ? x.value[1] : x.value

@doc doc"""
    RnTVector <: TVector

the tangent vector $\xi \in T_x\mathcal M$ for $\mathcal M=\mathbb R^n$
represented by an $n$-dimensional `Vector{T}`, where `T <: AbstractFloat`.
"""
struct RnTVector{T <: AbstractFloat}  <: TVector
  value::Vector{T}
  RnTVector{T}(value::Vector{T})  where {T <: AbstractFloat}  = new(value)
  RnTVector{T}(value::T) where {T <: AbstractFloat}  = new([value])
end
RnTVector(value::T) where {T <: AbstractFloat} = RnTVector{T}(value)
RnTVector(value::Vector{T})  where {T <: AbstractFloat}  = RnTVector{T}(value)

getValue(??::RnTVector) = length(??.value)==1 ? ??.value[1] : ??.value

# Traits
# ---
# (a) Rn is a MatrixManifold
@traitimpl IsMatrixM{Euclidean}
@traitimpl IsMatrixP{RnPoint}
@traitimpl IsMatrixTV{RnTVector}

# Functions
# ---
@doc doc"""
    distance(M,x,y)

compute the Euclidean distance $\lVert x - y\rVert$
"""
function distance(M::Euclidean,x::RnPoint{T},y::RnPoint{T})::T where {T <: AbstractFloat}
    if length(getValue(x)) > 1
        return norm( getValue(x) - getValue(y) )
    else
        return abs( getValue(x) - getValue(y) )
    end
end
@doc doc"""
    dot(M,x,??,??)
Computes the Euclidean inner product of `??` and `??`, i.e.
$\langle\xi,\nu\rangle = \displaystyle\sum_{k=1}^n \xi_k\nu_k$.
"""
dot(M::Euclidean,x::RnPoint{T},??::RnTVector{T}, ??::RnTVector{T}) where {T <: AbstractFloat} = dot( getValue(??) , getValue(??) )
@doc doc"""
    exp(M,x,??[, t=1.0])

compute the exponential map on the [`Euclidean`](@ref) manifold `M`, i.e.
$x+t*\xi$, where the scaling parameter `t` is optional.
"""
exp(M::Euclidean,x::RnPoint{T},??::RnTVector{T},t::Float64=1.0) where {T <: AbstractFloat} = RnPoint(getValue(x) + t*getValue(??) )
@doc doc"""
    injectivityRadius(M)

return the injectivity radius of the [`Euclidean`](@ref) manifold `M`$=\mathbb R^n$.
"""
injectivityRadius(M::Euclidean) = sqrt(M.dimension)
@doc doc"""
    log(M,x,y)

computes the logarithmic map on the [`Euclidean`](@ref) manifold `M`, i.e. $y-x$.
"""
log(M::Euclidean,x::RnPoint{T},y::RnPoint{T})  where {T <: AbstractFloat} = RnTVector( getValue(y) - getValue(x) )
@doc doc"""
    manifoldDimension(x)

return the manifold dimension of the [`RnPoint`](@ref) `x`, i.e. $n$.
"""
manifoldDimension(x::RnPoint) = length( getValue(x) )
@doc doc"""
    manifoldDimension(M)

return the manifold dimension of the [`Euclidean`](@ref) manifold `M`,
i.e. the length of the vectors stored in `M.dimension`, i.e. $n$.
"""
manifoldDimension(M::Euclidean) = M.dimension
@doc doc"""
    parallelTransport(M,x,y,??)

compute the parallel transport  the [`Euclidean`](@ref) manifold `M`, which is
the identity.
"""
parallelTransport(M::Euclidean, x::RnPoint{T}, y::RnPoint{T}, ??::RnTVector{T})  where {T <: AbstractFloat} = ??
@doc doc"""
    project(M,x,v)

project a $n$-dimensional `Vector{T}` v on the tangent space of the
[`RnPoint{T}`](@ref) `x`. Since the tangent space is identical to the
$\mathbb R^n$, the mapping can be realized with the identity, i.e.
$\operatorname{project}(M,x,v) = v$.
"""
project(M::Euclidean, x::RnPoint{T}, v::Vector{T}) where {T <: AbstractFloat} = RnTVector{T}(v)
@doc doc"""
    randomMPoint(M[,T=Float64])

generate a random point on the [`Euclidean`](@ref) manifold `M`, where the
optional parameter determines the type of the entries of the
resulting [`RnPoint`](@ref).
"""
randomMPoint(M::Euclidean, T::DataType=Float64) = RnPoint( randn(T,M.dimension) )

doc"""
    randomTVector(M,x,:Gaussian[,??=1.0])

generate a Gaussian random vector on the [`Euclidean`](@ref) manifold `M` with
standard deviation `??`.
"""
randomTVector(M::Euclidean, x::RnPoint{T}, ::Val{:Gaussian}, ??::Float64=1.0) where {T} = RnTVector( ?? * randn(T,M.dimension) )

@doc doc"""
    (??,??) = tangentONB(M,x,y)

compute an ONB within the tangent space $T_x\mathcal M$ at the [`MPoint`](@ref)
on the [`Euclidean`](@ref) manifold `M`, such that $\xi=\log_xy$
is the first vector and compute the eigenvalues of the curvature tensor
$R(\Xi,\dot g)\dot g$, where $g=g_{x,\xi}$ is the geodesic with $g(0)=x$,
$\dot g(0) = \xi$, i.e. $\kappa_1$ corresponding to $\Xi_1=\xi$ is zero.

# See also
 [`jacobiField`](@ref), [`adjointJacobiField`](@ref).
"""
tangentONB(M::Euclidean, x::RnPoint{T}, y::RnPoint{T}) where {T <: AbstractFloat}  = tangentONB(M,x,log(M,x,y))
tangentONB(M::Euclidean, x::RnPoint{T}, ??::RnTVector{T}) where {T <: AbstractFloat}  =
  [ RnTVector( Matrix{T}(I,manifoldDimension(x),manifoldDimension(x))[:,i] )
    for i in 1:manifoldDimension(x)], zeros(manifoldDimension(x))

typeofTVector(::Type{RnPoint{T}}) where T = RnTVector{T}
typeofMPoint(::Type{RnTVector{T}}) where T = RnPoint{T}

@doc doc"""
    typicalDistance(M)

returns the typical distance on the [`Euclidean`](@ref) manifold `M`: $\sqrt{n}$.
"""
typicalDistance(M::Euclidean) = sqrt(M.dimension)

@doc doc"""
    validateMPoint(M,x)

Checks that a [`RnPoint`](@ref) `x` has a valid value for a point on the
[`Euclidean`](@ref) manifold `M`$=\mathbb R^n$, which is the case if the
dimensions fit.
"""
validateMPoint(M::Euclidean, x::RnPoint) = manifoldDimension(M) == manifoldDimension(x)

@doc doc"""
    validateTVector(M,x,??)

Checks, that the [`RnTVector`](@ref) `??` is a valid tangent vector in the
tangent space of the [`RnPoint`](@ref) `x` ont the [`Euclidean`](@ref)
manifold `M`, which is always the case as long as their vector dimensions agree.
"""
validateTVector(M::Euclidean,x::RnPoint,??::RnTVector) = length(getValue(x) )== length(getValue(??))

@doc doc"""
    ?? = zeroTVector(M,x)

returns a zero vector in the tangent space $T_x\mathcal M$ of the
[`RnPoint`](@ref) $x\in\mathbb R^n$ on the [`Euclidean`](@ref) manifold `M`.
"""
function zeroTVector(M::Euclidean, x::RnPoint{T}) where {T <: AbstractFloat}
    return RnTVector(  zero( getValue(x) )  )
end
#
#
# --- Display functions for the objects/types
show(io::IO, M::Euclidean) = print(io, "The $(M.name)");
show(io::IO, x::RnPoint) = print(io, "Rn($( getValue(x) ))");
show(io::IO, ??::RnTVector) = print(io, "RnT($( getValue(??) ))");
