import Adapt

@inline fractional_index(x, xs, Δx) = @inbounds (x - xs[1]) / Δx

# Linear Lagrange polynomials
@inline ϕ₀(ξ) = 1 - ξ
@inline ϕ₁(ξ) = ξ

#####
##### Linear interpolation for regularly spaced vector data e.g. T(z) or τ(t).
#####

struct InterpolatedTimeSeries{D, T, Δ, M}
     data :: D
        t :: T
       Δt :: Δ
    t_min :: M
    t_max :: M
end

@inline function _interpolate(profile::InterpolatedTimeSeries, t★)
    n = 1 + fractional_index(t★, profile.t, profile.Δt)
    τ, n = mod(n, 1), Base.unsafe_trunc(Int, n)

    # Ensure we don't go out of bounds.
    if t★ <= profile.t_min
        return @inbounds profile.data[1]
    elseif t★ >= profile.t_max
        return @inbounds profile.data[n]
    else
        return @inbounds ϕ₀(τ) * profile.data[n] + ϕ₁(τ) * profile.data[n+1]
    end
end

@inline (profile::InterpolatedTimeSeries)(t) = _interpolate(profile, t)

@inline function ∂t(profile::InterpolatedTimeSeries, t★)
    n = 1 + fractional_index(t★, profile.t, profile.Δt)
    n = Base.unsafe_trunc(Int, n)

    # Ensure we don't go out of bounds.
    if t★ <= profile.t_min || t★ >= profile.t_max
        return 0
    else
        return @inbounds (profile.data[n+1] - profile.data[n]) / profile.Δt
    end
end

Adapt.adapt_structure(to, profile::InterpolatedTimeSeries) =
    InterpolatedTimeSeries(Adapt.adapt(to, profile.data),
                           Adapt.adapt(to, profile.t),
                           Adapt.adapt(to, profile.Δt),
                           Adapt.adapt(to, profile.t_min),
                           Adapt.adapt(to, profile.t_max))

#####
##### Linear interpolated for regularly spaced profile time series, e.g. T(z, t).
#####

struct InterpolatedProfileTimeSeries{D, Z, T, M}
     data :: D
        z :: Z
        t :: T
       Δz :: M
       Δt :: M
    z_max :: M
    t_max :: M
end

@inline function _interpolate(profile::InterpolatedProfileTimeSeries, z★, t★)
    k = 1 + fractional_index(z★, profile.z, profile.Δz)
    n = 1 + fractional_index(t★, profile.t, profile.Δt)

    ζ, k = mod(k, 1), Base.unsafe_trunc(Int, k)
    τ, n = mod(n, 1), Base.unsafe_trunc(Int, n)

    # Ensure we don't go out of bounds.
    if z★ >= profile.z_max && t★ >= profile.t_max
        return @inbounds profile.data[k, n]
    elseif z★ >= profile.z_max
        return @inbounds ϕ₀(τ) * profile.data[k, n] + ϕ₁(τ) * profile.data[k, n+1]
    elseif t★ >= profile.t_max
        return @inbounds ϕ₀(ζ) * profile.data[k, n] + ϕ₁(ζ) * profile.data[k+1, n]
    else
        return @inbounds (  ϕ₀(ζ) * ϕ₀(τ) * profile.data[k,   n  ]
                          + ϕ₀(ζ) * ϕ₁(τ) * profile.data[k,   n+1]
                          + ϕ₁(ζ) * ϕ₀(τ) * profile.data[k+1, n  ]
                          + ϕ₁(ζ) * ϕ₁(τ) * profile.data[k+1, n+1])
    end
end

@inline (profile::InterpolatedProfileTimeSeries)(z, t) = _interpolate(profile, z, t)

@inline function ∂z(profile::InterpolatedProfileTimeSeries, z★, t★)
    k = 1 + fractional_index(z★, profile.z, profile.Δz)
    n = 1 + fractional_index(t★, profile.t, profile.Δt)

    ζ, k = mod(k, 1), Base.unsafe_trunc(Int, k)
    τ, n = mod(n, 1), Base.unsafe_trunc(Int, n)

    # Ensure we don't go out of bounds.
    z★ >= profile.z_max && (k = k - 1)
    t★ >= profile.t_max && (n = n - 1)

    dataᵏ   = @inbounds ϕ₀(τ) * profile.data[k,   n] + ϕ₁(τ) * profile.data[k,   n+1]
    dataᵏ⁺¹ = @inbounds ϕ₀(τ) * profile.data[k+1, n] + ϕ₁(τ) * profile.data[k+1, n+1]

    return (dataᵏ⁺¹ - dataᵏ) / profile.Δz
end

@inline function ∂t(profile::InterpolatedProfileTimeSeries, z★, t★)
    k = 1 + fractional_index(z★, profile.z, profile.Δz)
    n = 1 + fractional_index(t★, profile.t, profile.Δt)

    ζ, k = mod(k, 1), Base.unsafe_trunc(Int, k)
    τ, n = mod(n, 1), Base.unsafe_trunc(Int, n)

    # Ensure we don't go out of bounds.
    z★ >= profile.z_max && (k = k - 1)
    t★ >= profile.t_max && (n = n - 1)

    dataⁿ   = @inbounds ϕ₀(ζ) * profile.data[k,   n] + ϕ₁(ζ) * profile.data[k+1,   n]
    dataⁿ⁺¹ = @inbounds ϕ₀(ζ) * profile.data[k, n+1] + ϕ₁(ζ) * profile.data[k+1, n+1]

    return (dataⁿ⁺¹ - dataⁿ) / profile.Δt
end

Adapt.adapt_structure(to, profile::InterpolatedProfileTimeSeries) =
    InterpolatedProfileTimeSeries(Adapt.adapt(to, profile.data),
                                  Adapt.adapt(to, profile.z),
                                  Adapt.adapt(to, profile.t),
                                  Adapt.adapt(to, profile.Δz),
                                  Adapt.adapt(to, profile.Δt),
                                  Adapt.adapt(to, profile.z_max),
                                  Adapt.adapt(to, profile.t_max))
