using Oceananigans

using Oceananigans.Utils: prettytime

function experiment_name(config; prefix, postfix="", overwrite_existing=false)
    s = prefix

    s *= "_" * replace(prettytime(config[:stop_time]), " " => "")
    s *= postfix

    if isdir(s) && !overwrite_existing
        n = 2
        while isdir(s * string(n))
            n = n + 1
        end
        s *= string(n)
    end

    mkpath(s)

    return s
end

config = Dict(

    :directories => Dict(
        :ecco_data_dir => "/storage6/alir/ecco4v4",
        :figures_dir => "/storage6/alir/RealityInspiredOceanTurbulary/site_South_of_Durban_gpu_test4/figures",
        :simulation_output_dir => "/storage6/alir/RealityInspiredOceanTurbulary/site_South_of_Durban_gpu_test4/data",
    ),

    :constants => Dict(
        :rotation_rate => 7.2921e-5,
        :gravitational_acceleration => 9.81,

        # Reference density of seawater [kg/m³] used by ECCO's equation of state (UNESCO).
        :reference_density_ecco => 1029,

        # Specific heat capacity of seawater at constant pressure [J/(kg·K)].
        :specific_heat_capacity_seawater => 4000
    ),

    :site => Dict(
        :name => "South_of_Durban",
        :latitude => -35,
        :longitude => 32,
        :start_date => Date(2014, 01, 01),
        :end_date => Date(2014, 12, 31)
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
