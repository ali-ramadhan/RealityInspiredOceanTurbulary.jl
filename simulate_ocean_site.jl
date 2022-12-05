using Logging
using Printf
using Dates

using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials.TEOS10

using Oceananigans.Architectures: array_type
using Oceananigans.BuoyancyModels: BuoyancyField

config_file = ARGS[1]
include(config_file)

Logging.global_logger(OceananigansLogger())

include("ocean_site_ecco_data.jl")
include("ocean_site_analysis.jl")
include("interpolated_profiles.jl")


@info "Gathering site data..."
site = gather_site_data()

if config[:plot_site_forcing]
    @info "Plotting site surface forcing..."
    plot_site_forcing(site)
end

if config[:animate_site_profiles]
    @info "Animating site profiles..."
    animate_site_profiles(site)
end


@info "Mapping grid..."

architecture = config[:architecture]
ArrayType = array_type(architecture)

topology = (Periodic, Periodic, Bounded)

Nx = config[:grid][:Nx]
Ny = config[:grid][:Ny]
Nz = config[:grid][:Nz]

Lx = config[:grid][:Lx]
Ly = config[:grid][:Ly]
Lz = config[:grid][:Lz]

domain = (x=(0, Lx), y=(0, Ly), z=(-Lz, 0))

grid = RectilinearGrid(architecture; topology, size=(Nx, Ny, Nz), domain...)

@info "Grid: $grid."


@info "Interpolating ECCO data..."

ℑτx = interpolate_ecco_timeseries(site, site["EXFtaue"]; ArrayType)
ℑτy = interpolate_ecco_timeseries(site, site["EXFtaun"]; ArrayType)
ℑQΘ = interpolate_ecco_timeseries(site, site["TFLUX"]; ArrayType)
ℑQS = interpolate_ecco_timeseries(site, site["SFLUX"]; ArrayType)

regular_zgrid = RectilinearGrid(architecture; topology=(Flat, Flat, Bounded), size=Nz, z=(-Lz, 0))
ecco_zgrid = ecco_vertical_grid()

ℑU = interpolate_profile(site["EVEL"], site["time"], ecco_zgrid, regular_zgrid; ArrayType)
ℑV = interpolate_profile(site["NVEL"], site["time"], ecco_zgrid, regular_zgrid; ArrayType)
ℑΘ = interpolate_profile(site["THETA"], site["time"], ecco_zgrid, regular_zgrid; ArrayType)
ℑS = interpolate_profile(site["SALT"], site["time"], ecco_zgrid, regular_zgrid; ArrayType)
ℑU_geo = interpolate_profile(site["U_geo"], site["time"], ecco_zgrid, regular_zgrid; ArrayType)
ℑV_geo = interpolate_profile(site["V_geo"], site["time"], ecco_zgrid, regular_zgrid; ArrayType)

# fig = Figure()
# ax = Axis(fig[1, 1])
# lines!(ax, interior(f_ecco)[:], znodes(Center, ecco_zgrid))
# scatter!(ax, ℑU.(znodes(Center, regular_zgrid), 0), znodes(Center, regular_zgrid), marker='o')
# ylims!(ax, (-512, 0))
# save("tmp.png", fig, px_per_unit=3)
