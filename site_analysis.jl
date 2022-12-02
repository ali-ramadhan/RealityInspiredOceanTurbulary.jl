using TOML
using Dates
using Match
using ProgressMeter
using NCDatasets

config = TOML.parsefile("configuration.toml")

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

function ecco_dataset(var, date)
    datestr = Dates.format(date, "YYYY-mm-dd")

    filename = @match var begin
        "MXLDEPTH" => "OCEAN_MIXED_LAYER_DEPTH_day_mean_$(datestr)_ECCO_V4r4_latlon_0p50deg.nc"
        "TFLUX" => "OCEAN_AND_ICE_SURFACE_HEAT_FLUX_day_mean_$(datestr)_ECCO_V4r4_latlon_0p50deg.nc"
        "SFLUX" => "OCEAN_AND_ICE_SURFACE_FW_FLUX_day_mean_$(datestr)_ECCO_V4r4_latlon_0p50deg.nc"
        "EVEL" || "NVEL" || "WVEL" => "OCEAN_VELOCITY_day_mean_$(datestr)_ECCO_V4r4_latlon_0p50deg.nc"
        "THETA" || "SALT" => "OCEAN_TEMPERATURE_SALINITY_day_mean_$(datestr)_ECCO_V4r4_latlon_0p50deg.nc"
    end

    filepath = joinpath(config["directories"]["ecco_data_dir"], filename)
    return NCDataset(filepath)
end

function _ecco_state(ds, var, idx_lat, idx_lon)
    return @match var begin
        "MXLDEPTH" || "TFLUX" || "SFLUX" => ds[var][idx_lon, idx_lat, 1]
        "EVEL" || "NVEL" || "WVEL" || "THETA" || "SALT" => ds[var][idx_lon, idx_lat, :, 1]
    end
end

function ecco_state(var, lat, lon, date)
    ds = ecco_dataset(var, date)
    idx_lat, idx_lon = ecco_indices(ds, lat, lon)
    return _ecco_state(ds, var, idx_lat, idx_lon)
end

function ecco_state_timeseries(var, lat, lon, start_date, end_date)
    date_range = start_date:Day(1):end_date
    return @showprogress 0.1 "$var..." [ecco_state(var, lat, lon, date) for date in date_range]
end

function concat_vectors(ts)
    if ts[1] isa Vector
        return hcat(ts...)
    else
        return ts
    end
end

site_lat = config["site"]["latitude"]
site_lon = config["site"]["longitude"]
start_date = config["site"]["start"]
end_date = config["site"]["end"]

@info "Getting ECCO state @ ($(site_lat)°N, $(site_lon)°E) from $start_date to $end_date..."

ecco_vars = ["MXLDEPTH", "TFLUX", "SFLUX", "EVEL", "NVEL", "WVEL", "THETA", "SALT"]

site = Dict(
    var => concat_vectors(ecco_state_timeseries(var, site_lat, site_lon, start_date, end_date))
    for var in ecco_vars
)

site["time"] = collect(start_date:Day(1):end_date)


