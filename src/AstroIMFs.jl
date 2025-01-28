module AstroIMFs

using QuadGK
using Unitful, UnitfulAstro

export AbstractIMF,
    PowerlawIMF,
    eval_by_mass,
    eval_by_num,
    integrate_eval_by_num,
    integrate_eval_by_mass,
    mass_in_massinterval,
    nstars_in_massinterval

"""
    abstract type AbstractIMF end

The abstract type for stellar initial mass functions.

The following functions can be used on objects of this supertype:
- [`eval_by_mass`](@ref)
- [`eval_by_num`](@ref)
- [`integrate_eval_by_mass`](@ref)
- [`integrate_eval_by_num`](@ref)
- [`mass_in_massinterval`](@ref)
- [`nstars_in_massinterval`](@ref)

# For developers
The following methods have to be defined for subtypes:
- [`eval_by_mass`](@ref)
- [`eval_by_num`](@ref)
"""
abstract type AbstractIMF end
Broadcast.broadcastable(imf::AbstractIMF) = Ref(imf)


"""
    mass_in_massinterval(imf::AbstractIMF, mtot::Number, mrange)

Return the expected total mass found within the mass range `mrange` for the `imf` of mass `mtot`.

The mass range is given as a Tuple or AbstractArray.
"""
function mass_in_massinterval(imf::AbstractIMF, mtot::Number, mrange)
    int, = integrate_eval_by_mass(imf, mrange)
    return int * mtot
end


"""
    nstars_in_massinterval(imf::AbstractIMF, mtot::Number, mrange)

Return the expected total number of stars found within the mass range `mrange` for the `imf` of mass `mtot`.

The mass range is given as a Tuple or AbstractArray.
"""
function nstars_in_massinterval(imf::AbstractIMF, mtot::Number, mrange)
    int, = integrate_eval_by_num(imf, mrange)
    return int * mtot
end


"""
    struct PowerlawIMF{M,S,T} <: AbstractIMF
        massmin::M
        massmax::M
        slopemasses::Vector{M}
        slopes::Vector{S}
        slopenorms::Vector{T}
    end

A powerlaw IMF with multiple separate powerlaw sections ranging the mass range from `massmin` to `massmax`.

This IMF can also handle unitful quantities, i.e. using solar masses as units.

# Attributes
- `slopemasses`: vector of edge masses between the slopes of length `N`, where the last element is equal
  to `massmin` (and the first element must be smaller or equal to `massmax`)
- `slopes`: vector of the powerlaw segment slopes
- `slopenorms`: vector of normalization values for the powerlaw segments to force the integral to be 1


# Examples
```julia
# for a Chabrier powerlaw IMF:
imf_chabrier = PowerlawIMF(0.1, 100.0, [1, .3, .1], [1.3, 0.8, 0.2], [0.241367, 0.241367, 0.497056])
```
"""
struct PowerlawIMF{M,S,T} <: AbstractIMF
    massmin::M
    massmax::M
    slopemasses::Vector{M}
    slopes::Vector{S}
    slopenorms::Vector{T}

    function PowerlawIMF(
        massmin::M,
        massmax::M,
        slopemasses::Vector{M},
        slopes::Vector{S},
        slopenorms::Vector{T},
    ) where {M,S,T}
        @assert length(slopemasses) == length(slopes) == length(slopenorms)
        if length(slopes) > 1
            @assert slopemasses[1] ≤ massmax
            @assert slopemasses[end] == massmin
            @assert all(slopenorms .< 1)
        end
        imf = new{M,S,T}(massmin, massmax, slopemasses, slopes, slopenorms)

        # verify that the IMF is properly normalized
        @assert isapprox(integrate_eval_by_mass(imf, (massmin, massmax)), 1; atol=0.001)

        return imf
    end
end

_ustrip(m::Number) = m
_ustrip(m::Unitful.Mass) = ustrip(u"Msun", m)

"""
    eval_by_mass(imf::PowerlawIMF, m::Number)

Returns the IMF probability density function evaluated at the mass `m`.

Integrating this over all masses results in 1.
"""
function eval_by_mass(imf::PowerlawIMF{M,S,T}, m::Number) where {M,S,T}
    imf.massmin ≤ m ≤ imf.massmax || return zero(T)

    i = searchsortedfirst(imf.slopemasses, m; rev=true)
    return imf.slopenorms[i] * _ustrip(m)^-imf.slopes[i]
end
eval_by_mass(imf::PowerlawIMF, m::Unitful.Mass) = @invoke(eval_by_mass(imf::PowerlawIMF, m::Number)) / u"Msun"

"""
    eval_by_num(imf::PowerlawIMF, m::Number)

Returns the IMF probability density function divided by the mass evaluated at the mass `m`.

Multiplying the integral of this over a mass interval with the total mass of stars that the IMF applies to
results in the expected number of stars within that interval.
"""
function eval_by_num(imf::PowerlawIMF{M,S,T}, m::Number) where {M,S,T}
    imf.massmin ≤ m ≤ imf.massmax || return zero(T)

    i = searchsortedfirst(imf.slopemasses, m; rev=true)
    return imf.slopenorms[i] * _ustrip(m)^-(1 + imf.slopes[i])
end
eval_by_num(imf::PowerlawIMF, m::Unitful.Mass) = @invoke(eval_by_num(imf::PowerlawIMF, m::Number)) / u"Msun^2"


function _partial_integrated_eval_by_num(slopenorm, slope, inf, sup)
    slopenorm * (1 / slope * (inf^-slope - sup^-slope))
end

function _partial_integrated_eval_by_mass(slopenorm, slope, inf, sup)
    slopenorm * (1 / (1 - slope) * (sup^(1 - slope) - inf^(1 - slope)))
end

"""
    integrate_eval_by_num(imf::PowerlawIMF, mrange)

Return the integral over [`eval_by_num`](@ref) for the mass range `mrange` as a Tuple or AbstractArray.
"""
function integrate_eval_by_num(imf::AbstractIMF, mrange)
    func = let imf = imf
        m -> eval_by_num(imf, m)
    end

    return quadgk(func, mrange[1], mrange[2])
end

"""
    integrate_eval_by_mass(imf::PowerlawIMF, mrange)

Return the integral over [`eval_by_mass`](@ref) for the mass range `mrange` as a Tuple or AbstractArray.
"""
function integrate_eval_by_mass(imf::AbstractIMF, mrange)
    func = let imf = imf
        m -> eval_by_mass(imf, m)
    end

    return quadgk(func, mrange[1], mrange[2])
end

function integrate_eval_by_num(imf::PowerlawIMF, mrange)
    _integrate_eval_by_partial(imf, mrange, _partial_integrated_eval_by_num)
end

function integrate_eval_by_num(imf::PowerlawIMF{M}, mrange) where {M<:Unitful.Mass}
    _integrate_eval_by_partial(imf, mrange, _partial_integrated_eval_by_num) / u"Msun"
end


"""
    integrate_eval_by_num(imf::PowerlawIMF, mrange)

Return the integral over [`eval_by_mass`](@ref) for the mass range `mrange` as a Tuple or AbstractArray.
"""
function integrate_eval_by_mass(imf::PowerlawIMF, mrange)
    _integrate_eval_by_partial(imf, mrange, _partial_integrated_eval_by_mass)
end

function _integrate_eval_by_partial(imf::PowerlawIMF{M,S,T}, mrange, partial_func) where {M,S,T}
    mrange = clamp.(mrange, imf.massmin, imf.massmax)
    mrange[1] == mrange[2] && return zero(T)

    int = zero(T)
    massmaxi = imf.massmax
    for i in eachindex(imf.slopes)
        (mrange[1] ≥ massmaxi || mrange[2] ≤ imf.slopemasses[i]) && continue # no overlap

        slope = imf.slopes[i]
        inf = max(mrange[1], imf.slopemasses[i])
        sup = min(mrange[2], massmaxi)
        int += partial_func(imf.slopenorms[i], slope, _ustrip(inf), _ustrip(sup))

        massmaxi = imf.slopemasses[i]
    end
    return int
end

end
