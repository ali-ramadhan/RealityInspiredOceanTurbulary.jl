config = Dict(

    :directories => Dict(
        :ecco_data_dir => "/storage6/alir/ecco4v4",
        :figures_dir => "/storage6/alir/RealityInspiredOceanTurbulary"
    ),

    :constants => Dict(
        :gravitational_acceleration => 9.81,
        :ecco_reference_density => 1029,
        :rotation_rate => 7.2921e-5
    ),

    :site => Dict(
        :name => "South of Durban",
        :latitude => -50,
        :longitude => 30,
        :start_date => Date(2014, 01, 01),
        :end_date => Date(2014, 12, 31)
    )

)
