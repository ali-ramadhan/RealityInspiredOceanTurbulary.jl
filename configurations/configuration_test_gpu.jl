using Oceananigans

using Dates: format
using Oceananigans.Utils: prettytime

config = Dict(

    :directories => Dict(
        :ecco_data_dir => "/storage6/alir/ecco4v4"
    ),

    :constants => Dict(
        :rotation_rate => 7.2921e-5,
        :gravitational_acceleration => 9.81,

        # Reference density of seawater [kg/m³] used by ECCO's equation of state (UNESCO).
        :reference_density_ecco => 1029,

        # Specific heat capacity of seawater at constant pressure [J/(kg·K)].
        :specific_heat_capacity_seawater => 4000,

        # Reference salinity [psu] used to compute the salt flux from E - P - R.
        :reference_salinity = 40
    ),

    :site => Dict(
        :name => "South_of_Durban",
        :latitude => -35,
        :longitude => 32,
        :start_date => Date(2014, 01, 01),
        :end_date => Date(2014, 01, 31)
    ),

    :architecture => GPU(),

    :grid => Dict(
        :Nx => 128,
        :Ny => 128,
        :Nz => 128,
        :Lx => 500,
        :Ly => 500,
        :Lz => 500
    ),

    :plot_site_forcing => true,
    :animate_site_profiles => true
)

# Determine ouptut directories

__BASE_OUTPUT_DIR = "/storage6/alir/RealityInspiredOceanTurbulary"

function experiment_dir(config; prefix="", postfix="", base_dir=__BASE_OUTPUT_DIR, overwrite_existing=false)
    s = prefix

    site_name = config[:site][:name]
    site_lat = config[:site][:latitude]
    site_lon = config[:site][:longitude]
    site_start = format(config[:site][:start_date], "YYYY-mm-dd")
    site_end = format(config[:site][:end_date], "YYYY-mm-dd")

    s *= "$(site_name)_lat$(site_lat)_lon$(site_lon)_$(site_start)_to_$(site_end)"

    s *= postfix

    dir = joinpath(base_dir, s)

    if isdir(dir) && !overwrite_existing
        n = 2
        while isdir(dir * string(n))
            n = n + 1
        end
        dir *= string(n)
    end

    mkpath(dir)

    return dir
end

__output_dir = joinpath(__BASE_OUTPUT_DIR, experiment_dir(config, postfix="_test_gpu"))

config[:directories][:figures_dir] = joinpath(__output_dir, "figures")
config[:directories][:simulation_output_dir] = joinpath(__output_dir, "data")
