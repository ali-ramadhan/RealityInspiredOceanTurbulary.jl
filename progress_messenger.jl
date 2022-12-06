using Printf
using Statistics

using Dates: format

mutable struct SimulationProgressMessenger{T, F} <: Function
      wall_time₀ :: T  # Wall time at simulation start
      wall_time⁻ :: T  # Wall time at previous calback
    compile_time :: T  # Estimate of simulation compile time
    cfl_tapering :: F  # Function to modify max CFL as a function of time
end

SimulationProgressMessenger(; cfl_tapering=nothing) =
SimulationProgressMessenger(
                  1e-9 * time_ns(),
                  1e-9 * time_ns(),
                  0.0,
                  cfl_tapering)

function (pm::SimulationProgressMessenger)(simulation)
    model = simulation.model
    i, t = model.clock.iteration, model.clock.time
    iteration_interval = simulation.callbacks[:progress].schedule.interval

    if !isnothing(pm.cfl_tapering)
        wizard = simulation.callbacks[:wizard].func
        wizard.cfl = cfl_tapering(t)
        wizard.diffusive_cfl = cfl_tapering(t)
    end

    # For more accurate ETA calculations, assume the first iteration interval
    # took the same amount of time as the second, so the remainder must be
    # compile time.
    if i == 2iteration_interval
        runtime_first = pm.wall_time⁻ - pm.wall_time₀
        runtime_second = 1e-9 * time_ns() - pm.wall_time⁻
        if runtime_first > runtime_second
            pm.compile_time = runtime_first - runtime_second
        end
    end

    progress = t / simulation.stop_time
    current_wall_time = 1e-9 * time_ns() - pm.wall_time₀

    time_since_last_callback = 1e-9 * time_ns() - pm.wall_time⁻
    wall_time_per_step = time_since_last_callback / iteration_interval
    pm.wall_time⁻ = 1e-9 * time_ns()

    u_max = maximum(abs, model.velocities.u)
    v_max = maximum(abs, model.velocities.v)
    w_max = maximum(abs, model.velocities.w)

    T_min = minimum(model.tracers.T)
    T_max = maximum(model.tracers.T)
    T_mean = mean(model.tracers.T)

    cfl_adv = AdvectiveCFL(simulation.Δt)(simulation.model)
    cfl_dif = DiffusiveCFL(simulation.Δt)(simulation.model)

    msg_line1 = @sprintf("[%06.2f%%] iteration: % 6d, time: % 10s, Δt: % 10s, wall time: % 8s (% 8s / time step)",
                          100 * progress, i, prettytime(t),
                          prettytime(simulation.Δt),
                          prettytime(current_wall_time),
                          prettytime(wall_time_per_step))

    # ETA calculation only accurate after 2 iteration intervals.
    if i >= 2iteration_interval
        ETA = (1 - progress) / progress * (current_wall_time - pm.compile_time)
        ETA_datetime = now() + Second(round(Int, ETA))
        eta_msg = @sprintf(", ETA: %s (%s)", format(ETA_datetime, "yyyy-mm-dd HH:MM:SS"), prettytime(ETA))
        msg_line1 *= eta_msg
    end

    @info msg_line1

    @info @sprintf("          └── CFL: %.2e, νκCFL: %.2e,  u⃗_max: (%.2e, %.2e, %.2e) m/s, T: (min=%.2f, mean=%.2f, max=%.2f)",
                    cfl_adv, cfl_dif, u_max, v_max, w_max, T_min, T_mean, T_max)

    @info ""

    return nothing
end
