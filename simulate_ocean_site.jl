using Logging
using Printf
using Dates

using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials.TEOS10

using CUDA: CuArray
using Oceananigans.Architectures: array_type
using Oceananigans.BuoyancyModels: BuoyancyField

config_file = ARGS[1]
include(config_file)

Logging.global_logger(OceananigansLogger())

include("ocean_site_ecco_data.jl")
include("ocean_site_analysis.jl")
include("interpolated_profiles.jl")
include("interpolate_ecco_data.jl")
include("progress_messenger.jl")
include("mixed_layer_depth.jl")
include("plot_les_solution.jl")


@info "Summoning ECCO data and diagnosing geostrophic background state..."

site = gather_site_data()

if config[:plot_site_forcing]
    @info "Plotting site surface forcing..."
    plot_site_forcing(site)
end

if config[:animate_site_profiles]
    @info "Animating site profiles..."
    animate_site_profiles(site)
end


@info "Selecting an architecture..."

architecture = config[:architecture]

ArrayType = architecture isa GPU ? CuArray{Float64} : Array{Float64}


@info "Mapping grid..."

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

# We use interpolated profiles on the CPU to set! field initial conditions and output to disk.

ℑτx_cpu = interpolate_ecco_timeseries(site, site["EXFtaue"], ArrayType=Array{Float64})
ℑτy_cpu = interpolate_ecco_timeseries(site, site["EXFtaun"], ArrayType=Array{Float64})
ℑQΘ_cpu = interpolate_ecco_timeseries(site, site["TFLUX"], ArrayType=Array{Float64})
ℑQS_cpu = interpolate_ecco_timeseries(site, site["SFLUX"], ArrayType=Array{Float64})
ℑmld_cpu = interpolate_ecco_timeseries(site, site["MXLDEPTH"], ArrayType=Array{Float64})

ℑU_cpu = interpolate_profile(site["EVEL"], site["time"], ecco_zgrid, regular_zgrid, ArrayType=Array{Float64})
ℑV_cpu = interpolate_profile(site["NVEL"], site["time"], ecco_zgrid, regular_zgrid, ArrayType=Array{Float64})
ℑΘ_cpu = interpolate_profile(site["THETA"], site["time"], ecco_zgrid, regular_zgrid, ArrayType=Array{Float64})
ℑS_cpu = interpolate_profile(site["SALT"], site["time"], ecco_zgrid, regular_zgrid, ArrayType=Array{Float64})
ℑU_geo_cpu = interpolate_profile(site["U_geo"], site["time"], ecco_zgrid, regular_zgrid, ArrayType=Array{Float64})
ℑV_geo_cpu = interpolate_profile(site["V_geo"], site["time"], ecco_zgrid, regular_zgrid, ArrayType=Array{Float64})

@info "Forcing mean-flow interactions..."

@inline _U_geo(x, y, z, t, p) = p.ℑU_geo(z, t)
@inline _V_geo(x, y, z, t, p) = p.ℑV_geo(z, t)

U_geo = BackgroundField(_U_geo, parameters=(; ℑU_geo))
V_geo = BackgroundField(_V_geo, parameters=(; ℑV_geo))

background_fields = (u=U_geo, v=V_geo)


@info "Enforcing boundary conditions..."

## Set up boundary conditions to
##   1. impose wind stresses at the ocean surface, and
##   2. impose heat and salt fluxes at the ocean surface.

ρ₀ = config[:constants][:reference_density_ecco]
cₚ = config[:constants][:specific_heat_capacity_seawater]

@inline    wind_stress_x(x, y, t, p) =   p.ℑτx(t) / p.ρ₀
@inline    wind_stress_y(x, y, t, p) =   p.ℑτy(t) / p.ρ₀
@inline temperature_flux(x, y, t, p) = - p.ℑQΘ(t) / p.ρ₀ / p.cₚ
@inline        salt_flux(x, y, t, p) =   p.ℑQS(t) / p.ρ₀

u_wind_stress_bc = FluxBoundaryCondition(wind_stress_x, parameters=(; ℑτx, ρ₀))
v_wind_stress_bc = FluxBoundaryCondition(wind_stress_y, parameters=(; ℑτy, ρ₀))
T_surface_flux_bc = FluxBoundaryCondition(temperature_flux, parameters=(; ℑQΘ, ρ₀, cₚ))
S_surface_flux_bc = FluxBoundaryCondition(salt_flux, parameters=(; ℑQS, ρ₀))

u_bcs = FieldBoundaryConditions(top=u_wind_stress_bc)
v_bcs = FieldBoundaryConditions(top=v_wind_stress_bc)
T_bcs = FieldBoundaryConditions(top=T_surface_flux_bc)
S_bcs = FieldBoundaryConditions(top=S_surface_flux_bc)

boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs)


@info "Framing the model..."

model = NonhydrostaticModel(
    grid = grid,
    # clock = Clock(time=DateTime(config[:site][:start_date])),
    timestepper = :RungeKutta3,
    advection = architecture isa CPU ? UpwindBiasedThirdOrder() : WENO(),
    background_fields = background_fields,
    buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState()),
    tracers = (:T, :S),
    coriolis = FPlane(latitude=site["latitude"]),
    closure = AnisotropicMinimumDissipation(),
    boundary_conditions = boundary_conditions
)

@info "Model: $model."


@info "Initializing site..."

ε(σ) = σ * randn() # noise

U₀(x, y, z) = 0
V₀(x, y, z) = 0
W₀(x, y, z) = ε(1e-10)
Θ₀(x, y, z) = ℑΘ_cpu(z, 0)
S₀(x, y, z) = ℑS_cpu(z, 0)

set!(model, u=U₀, v=V₀, w=W₀, T=Θ₀, S=S₀)


@info "Spinning up the simulation..."

simulation = Simulation(model, Δt=1seconds, stop_time=1days)

wizard = TimeStepWizard(cfl=0.5, diffusive_cfl=0.25, max_Δt=1minutes, min_Δt=1e-3*seconds)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

simulation.callbacks[:progress] = Callback(SimulationProgressMessenger(), IterationInterval(20))

@info "Simulation: $simulation"


@info "Setting up output writers..."

simulation_output_dir = config[:directories][:simulation_output_dir]
mkpath(simulation_output_dir)

b = BuoyancyField(model)
output_fields = merge(model.velocities, model.tracers, (; b))

simulation.output_writers[:fields] =
    JLD2OutputWriter(model, output_fields;
        dir = simulation_output_dir,
        schedule = TimeInterval(1hours),
        filename = "ocean_site_fields.jld2",
        with_halos = true,
        overwrite_existing = true
    )

simulation.output_writers[:surface] =
    JLD2OutputWriter(model, output_fields;
        dir = simulation_output_dir,
        filename = "ocean_site_surface.jld2",
        indices = (:, :, Nz),
        schedule = TimeInterval(10minutes),
        with_halos = true,
        overwrite_existing = true
    )

simulation.output_writers[:checkpointer] =
    Checkpointer(model;
        dir = simulation_output_dir,
        prefix = "model_checkpoint",
        schedule = TimeInterval(5days)
    )

profiles = (
    U = Average(model.velocities.u, dims=(1, 2)),
    V = Average(model.velocities.v, dims=(1, 2)),
    T = Average(model.tracers.T, dims=(1, 2)),
    S = Average(model.tracers.S, dims=(1, 2)),
    B = Average(b, dims=(1, 2))
)

simulation.output_writers[:profiles] =
    JLD2OutputWriter(model, profiles;
        dir = simulation_output_dir,
        filename = "ocean_site_profiles.jld2",
        schedule = TimeInterval(10minutes),
        with_halos = true,
        overwrite_existing = true
    )

large_scale_outputs = Dict(
       "τx" => model -> ℑτx_cpu(model.clock.time),
       "τy" => model -> ℑτy_cpu(model.clock.time),
       "QΘ" => model -> ℑQΘ_cpu(model.clock.time),
       "QS" => model -> ℑQS_cpu(model.clock.time),
        "U" => model -> ℑU_cpu.(znodes(Center, model.grid)[:], model.clock.time),
        "V" => model -> ℑV_cpu.(znodes(Center, model.grid)[:], model.clock.time),
        "T" => model -> ℑΘ_cpu.(znodes(Center, model.grid)[:], model.clock.time),
        "S" => model -> ℑS_cpu.(znodes(Center, model.grid)[:], model.clock.time),
    "U_geo" => model -> ℑU_geo_cpu.(znodes(Center, model.grid)[:], model.clock.time),
    "V_geo" => model -> ℑV_geo_cpu.(znodes(Center, model.grid)[:], model.clock.time),
    "mld_ECCO" => model -> ℑmld_cpu(model.clock.time),
    "mld_LES" => model -> mixed_layer_depth(model)
)

large_scale_dims = Dict(
       "τx" => (),
       "τy" => (),
       "QΘ" => (),
       "QS" => (),
        "U" => ("zC",),
        "V" => ("zC",),
        "T" => ("zC",),
        "S" => ("zC",),
    "U_geo" => ("zC",),
    "V_geo" => ("zC",),
    "mld_ECCO" => (),
    "mld_LES" => (),
)

simulation.output_writers[:large_scale_nc] =
    NetCDFOutputWriter(model, large_scale_outputs;
        dir = simulation_output_dir,
        filename = "ocean_site_large_scale.nc",
        schedule = TimeInterval(10minutes),
        dimensions = large_scale_dims
    )



wave = raw"""
           _.====.._
         ,:._       ~-_
             `\        ~-_
               |          `.
             ,/             ~-_
..__-..__..-''                 ~~--..__...----... LESbrary.jl ..__..--
"""

fish = raw"""
o                      o                    o
  o                     o                   o
 o                     o                    o
o   .''''.            o   .''''.            o   .''''.
 o /O)    './|         o /O)    './|         o /O)    './|
   > ) \| .'\|           > ) \| .'\|           > ) \| .'\|
    `....`                `....`                `....`
      ` `                   ` `                   ` `
             o                      o                   o
            o                      o                    o
            o   .''''.            o   .''''.             o  .''''.
             o /O)    './|         o /O)    './|         o /O)    './|
               > ) \| .'\|           > ) \| .'\|           > ) \| .'\|
                `....`                `....`                `....`
                  ` `                   ` `                   ` `
"""

print(wave)
print(fish)


@info "Teaching the simulation to run!..."

run!(simulation)

ds_profiles = FieldDataset(simulation.output_writers[:profiles].filepath)
animate_les_profiles(ds_profiles)
