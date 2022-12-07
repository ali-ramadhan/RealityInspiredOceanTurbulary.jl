using Dates
using Match
using ProgressMeter
using NCDatasets
using Oceananigans

config_file = ARGS[1]
include(config_file)

function validate_lat_lon(lat, lon)
    @assert -90 <= lat <= 90
    @assert -180 <= lon <= 180
    return nothing
end

function ecco_indices(ds, lat, lon)
    lats = ds["latitude"][:]
    lons = ds["longitude"][:]

    _, idx_lat = findmin(abs, lats .- lat)
    _, idx_lon = findmin(abs, lons .- lon)

    return idx_lat, idx_lon
end

function ecco_dataset(var, date; dir=config[:directories][:ecco_data_dir])
    datestr = Dates.format(date, "YYYY-mm-dd")

    filename = @match var begin
        "EXFtaue" || "EXFtaun" => "OCEAN_AND_ICE_SURFACE_STRESS_day_mean_$(datestr)_ECCO_V4r4_latlon_0p50deg.nc"
        "TFLUX" => "OCEAN_AND_ICE_SURFACE_HEAT_FLUX_day_mean_$(datestr)_ECCO_V4r4_latlon_0p50deg.nc"
        "SFLUX" => "OCEAN_AND_ICE_SURFACE_FW_FLUX_day_mean_$(datestr)_ECCO_V4r4_latlon_0p50deg.nc"
        "MXLDEPTH" => "OCEAN_MIXED_LAYER_DEPTH_day_mean_$(datestr)_ECCO_V4r4_latlon_0p50deg.nc"
        "EVEL" || "NVEL" || "WVEL" => "OCEAN_VELOCITY_day_mean_$(datestr)_ECCO_V4r4_latlon_0p50deg.nc"
        "THETA" || "SALT" => "OCEAN_TEMPERATURE_SALINITY_day_mean_$(datestr)_ECCO_V4r4_latlon_0p50deg.nc"
        "RHOAnoma" => "OCEAN_DENS_STRAT_PRESS_day_mean_$(datestr)_ECCO_V4r4_latlon_0p50deg.nc"
    end

    filepath = joinpath(dir, filename)
    return NCDataset(filepath)
end

function _ecco_state(ds, var, idx_lat, idx_lon)
    return @match var begin
        "EXFtaue" || "EXFtaun" || "TFLUX" || "SFLUX" || "MXLDEPTH" => ds[var][idx_lon, idx_lat, 1]

        # Reverse the vertical profile so increasing the index increases the vertical coordinate (Oceananigans convention).
        "EVEL" || "NVEL" || "WVEL" || "THETA" || "SALT" || "RHOAnoma" => reverse(ds[var][idx_lon, idx_lat, :, 1])
    end
end

function ecco_state(var, lat, lon, date)
    ds = ecco_dataset(var, date)
    idx_lat, idx_lon = ecco_indices(ds, lat, lon)
    state = _ecco_state(ds, var, idx_lat, idx_lon)
    close(ds)
    return state
end

function ecco_state_timeseries(var, lat, lon, start_date, end_date)
    date_range = start_date:Day(1):end_date
    return @showprogress 0.1 "$var" [ecco_state(var, lat, lon, date) for date in date_range]
end

# Convert Vector{Vector{T}} into Matrix{T} for vertical profiles from ECCO.
function concat_vectors(ts)
    if ts[1] isa Vector
        return hcat(ts...)
    else
        return ts
    end
end

function ecco_grid(ds)
    latitude = ds["latitude"][:]
    longitude = ds["longitude"][:]
    z_faces = ds["Z_bnds"][:] |> unique |> sort

    topology = (Periodic, Bounded, Bounded)
    size = (length(longitude), length(latitude), length(z_faces)-1)
    grid = LatitudeLongitudeGrid(Float32; topology, size, longitude=(-180, 180), latitude=(-90, 90), z=z_faces)

    return grid
end

function ecco_vertical_grid()
    ds = ecco_dataset("EVEL", Date(2013, 01, 01))

    topology = (Flat, Flat, Bounded)

    z_faces = ds["Z_bnds"][:] |> unique |> sort
    Nz = length(z_faces) - 1

    grid = RectilinearGrid(size=Nz, z=z_faces, topology=topology, halo=1)

    close(ds)

    return grid
end

function _geostrophic_base_state(site, date)
    g = config[:constants][:gravitational_acceleration]
    ρ₀ = config[:constants][:reference_density_ecco]
    Ω = config[:constants][:rotation_rate]

    lat = site["latitude"]
    lon = site["longitude"]

    ds = ecco_dataset("RHOAnoma", date)
    idx_lat, idx_lon = ecco_indices(ds, lat, lon)

    grid = ecco_grid(ds)

    ρ′ = reverse(ds["RHOAnoma"][:][:, :, :, 1], dims=3)
    ρ′ = replace(ρ′, missing => NaN)  # Oceananigans Fields cannot handle `missing` values.

    B = Field{Center, Center, Center}(grid)
    set!(B, @. -g * ρ′ / ρ₀)

    ∂B∂x = @at (Center, Center, Center) ∂x(B)
    ∂B∂y = @at (Center, Center, Center) ∂y(B)

    ∂B∂x_site = ∂B∂x[idx_lon, idx_lat, 1:grid.Nz]
    ∂B∂y_site = ∂B∂y[idx_lon, idx_lat, 1:grid.Nz]

    ∂B∂x_site = replace(∂B∂x_site, NaN => 0)
    ∂B∂y_site = replace(∂B∂y_site, NaN => 0)

    land_points = ismissing.(site["EVEL"][:, 1])

    Δz = grid.Δzᵃᵃᶜ[1:grid.Nz]
    ∫dz_∂B∂x = cumsum(∂B∂x_site .* Δz)
    ∫dz_∂B∂y = cumsum(∂B∂y_site .* Δz)

    f = 2Ω * sind(lat)
    U_geo = -1/f * ∫dz_∂B∂y
    V_geo =  1/f * ∫dz_∂B∂x

    U_geo[land_points] .= NaN
    V_geo[land_points] .= NaN

    U_geo = replace(U_geo, NaN => missing)
    V_geo = replace(V_geo, NaN => missing)

    close(ds)

    return U_geo, V_geo
end

function geostrophic_base_state(site)
    U_geo = similar(site["EVEL"])
    V_geo = similar(site["NVEL"])

    Nt = length(site["time"])
    prog = Progress(Nt, desc="Diagnosing geostrophic base state...", showspeed=true)
    for (n, date) in enumerate(site["time"])
        U_geo_day, V_geo_day = _geostrophic_base_state(site, date)
        U_geo[:, n] .= U_geo_day
        V_geo[:, n] .= V_geo_day
        next!(prog, showvalues=[(:date, date)])
    end

    return U_geo, V_geo
end

function gather_site_data(;
        site_lat = config[:site][:latitude],
        site_lon = config[:site][:longitude],
        start_date = config[:site][:start_date],
        end_date = config[:site][:end_date]
    )

    @info "Getting ECCO state @ ($(site_lat)°N, $(site_lon)°E) from $start_date to $end_date..."

    site = Dict()

    site["time"] = collect(start_date:Day(1):end_date)
    site["latitude"] = site_lat
    site["longitude"] = site_lon

    dsU1 = ecco_dataset("EVEL", start_date)
    site["z"] = reverse(dsU1["Z"][:])
    close(dsU1)

    ecco_vars = ["EXFtaue", "EXFtaun", "TFLUX", "SFLUX", "MXLDEPTH", "EVEL", "NVEL", "THETA", "SALT", "RHOAnoma"]

    for var in ecco_vars
        site[var] = concat_vectors(ecco_state_timeseries(var, site_lat, site_lon, start_date, end_date))
    end

    U_geo, V_geo = geostrophic_base_state(site)

    site["U_geo"] = U_geo
    site["V_geo"] = V_geo

    return site
end
