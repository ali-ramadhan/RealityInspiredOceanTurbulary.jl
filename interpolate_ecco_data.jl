import Interpolations
using Interpolations: Gridded, Linear, interpolate, extrapolate

function interpolate_ecco_timeseries(site, timeseries; ArrayType=Array{Float64})
    times = site["time"]
    ℑt = [Second(t - times[1]).value for t in times]
    Δt = ℑt[2] - ℑt[1]
    t_min, t_max = extrema(ℑt)

    return InterpolatedTimeSeries(ArrayType(timeseries), ℑt, Δt, t_min, t_max)
end

function interpolate_profile(profile, times, ecco_grid, regular_grid; array_type=Array{Float64})
    ℑt = [Second(t - times[1]).value for t in times] |> array_type
    Δt = ℑt[2] - ℑt[1]
    t_max = last(ℑt)
    z_max = ecco_grid.zᵃᵃᶜ[ecco_grid.Nz]

    ## First we interpolate from the stretched ECCO grid onto the regular Oceananigans.jl grid.

    # Coordinates needs to be in increasing order for Interpolations.jl.
    zc_ecco = znodes(Center, ecco_grid)

    ℑprofile = interpolate((zc_ecco, ℑt), profile, Gridded(Linear()))

    # We might go out of bounds if the Oceananigans grid is finer in resolution. In that case, grab the nearest value.
    ℑprofile = extrapolate(ℑprofile, Interpolations.Flat())

    zc_regular = znodes(Center, regular_grid)[:]
    interpolated_data = ℑprofile.(zc_regular, ℑt')

    ## Then we construct and return an InterpolatedProfileTimeSeries for fast linear interpolation in kernels.
    return InterpolatedProfileTimeSeries(array_type(interpolated_data), zc_regular, ℑt, regular_grid.Δzᵃᵃᶜ, Δt, z_max, t_max)
end

function interpolate_profiles(sose_profiles, sose_grid, grid, times; array_type=Array{Float64})
    return (
        U    = interpolate_profile(sose_profiles.U,    sose_grid, grid, times, array_type),
        V    = interpolate_profile(sose_profiles.V,    sose_grid, grid, times, array_type),
        Θ    = interpolate_profile(sose_profiles.Θ,    sose_grid, grid, times, array_type),
        S    = interpolate_profile(sose_profiles.S,    sose_grid, grid, times, array_type),
        Ugeo = interpolate_profile(sose_profiles.Ugeo, sose_grid, grid, times, array_type),
        Vgeo = interpolate_profile(sose_profiles.Vgeo, sose_grid, grid, times, array_type)
    )
end
