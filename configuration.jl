using Oceananigans

config = Dict(

    :directories => Dict(
        :ecco_data_dir => "/storage6/alir/ecco4v4",
        :figures_dir => "/storage6/alir/RealityInspiredOceanTurbulary/site_South_of_Durban/figures",
        :simulation_output_dir => "/storage6/alir/RealityInspiredOceanTurbulary/site_South_of_Durban/data",
    ),

    :constants => Dict(
        :gravitational_acceleration => 9.81,
        :ecco_reference_density => 1029,
        :rotation_rate => 7.2921e-5
    ),

    :site => Dict(
        :name => "South_of_Durban",
        :latitude => -35,
        :longitude => 32,
        :start_date => Date(2014, 01, 01),
        :end_date => Date(2014, 12, 31)
    ),

    :architecture => CPU(),

    :grid => Dict(
        :Nx => 128,
        :Ny => 128,
        :Nz => 128,
        :Lx => 500,
        :Ly => 500,
        :Lz => 500
    ),

    :plot_site_forcing => false,
    :animate_site_profiles => false

)
