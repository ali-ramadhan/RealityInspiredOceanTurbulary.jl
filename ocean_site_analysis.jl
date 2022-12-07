using Dates
using Match
using ProgressMeter
using NCDatasets
using Oceananigans
using CairoMakie

config_file = ARGS[1]
include(config_file)

function plot_site_forcing(site; resolution=(1920, 1080), filepath=joinpath(config[:directories][:figures_dir], "site_forcing.png"))
    fig = Figure(; resolution)

    lat = site["latitude"]
    lon = site["longitude"]
    title = "Site analysis @ ($(lat)°N, $(lon)°E)"
    # Label(fig[0, :], title, fontsize=30)

    days = 0:length(site["time"])-1
    ndays = length(days)

    ax1 = Axis(fig[1, 1], ylabel="Wind stress")
    lines!(ax1, days, site["EXFtaue"], linewidth=3, label="τx (east)")
    lines!(ax1, days, site["EXFtaun"], linewidth=3, label="τy (north)")
    xlims!(ax1, (0, ndays))
    hidexdecorations!(ax1, grid=false)

    ax2 = Axis(fig[2, 1], ylabel="Surface heat flux (W/m²)")
    lines!(ax2, days, site["TFLUX"], linewidth=3)
    xlims!(ax2, (0, ndays))
    hidexdecorations!(ax2, grid=false)

    ax3 = Axis(fig[3, 1], ylabel="E - P - R (m/s)")
    lines!(ax3, days, site["oceFWflx"], linewidth=3)
    xlims!(ax3, (0, ndays))
    hidexdecorations!(ax3, grid=false)

    ax4 = Axis(fig[4, 1], xlabel="Time", ylabel="Mixed layer depth (m)", yreversed=true)
    lines!(ax4, days, site["MXLDEPTH"], linewidth=3)
    xlims!(ax4, (0, ndays))

    ax4.xtickformat = ds -> [Dates.format(site["time"][1] + Day(d), "YYYY-mm-dd") for d in ds]
    ax4.xticklabelrotation = π/4
    ax4.xticklabelalign = (:right, :center)

    linkxaxes!(ax1, ax2, ax3, ax4)

    @info "Plotting $filepath..."
    mkpath(dirname(filepath))
    save(filepath, fig)

    return nothing
end

function animate_site_profiles(site; resolution=(1920, 1080), framerate=15, filepath=joinpath(config[:directories][:figures_dir], "site_profiles.mp4"))
    time = site["time"]
    frames = 1:length(time)

    frame = Observable(1)

    lat = site["latitude"]
    lon = site["longitude"]
    title = @lift "Site analysis @ ($(lat)°N, $(lon)°E): $(time[$frame])"

    z = site["z"]
    U = @lift site["EVEL"][:, $frame]
    V = @lift site["NVEL"][:, $frame]
    T = @lift site["THETA"][:, $frame]
    S = @lift site["SALT"][:, $frame]
    ρ′ = @lift site["RHOAnoma"][:, $frame]

    U_geo = @lift site["U_geo"][:, $frame]
    V_geo = @lift site["V_geo"][:, $frame]

    fig = Figure(; resolution)

    ax1 = Axis(fig[1, 1], xlabel="U (m/s)", ylabel="z (m)")
    lines!(ax1, U, z, linewidth=3)
    lines!(ax1, U_geo, z, linewidth=3, linestyle=:dot)
    xlims!(ax1, extrema(skipmissing(site["EVEL"])))

    ax2 = Axis(fig[1, 2], xlabel="V (m/s)", ylabel="z (m)")
    lines!(ax2, V, z, linewidth=3)
    lines!(ax2, V_geo, z, linewidth=3, linestyle=:dot)
    xlims!(ax2, extrema(skipmissing(site["NVEL"])))
    hideydecorations!(ax2, grid=false)

    ax3 = Axis(fig[1, 3], xlabel="T (°C)", ylabel="z (m)")
    lines!(ax3, T, z, linewidth=3)
    xlims!(ax3, extrema(skipmissing(site["THETA"])))
    hideydecorations!(ax3, grid=false)

    ax4 = Axis(fig[1, 4], xlabel="S (psu)", ylabel="z (m)")
    lines!(ax4, S, z, linewidth=3)
    xlims!(ax4, extrema(skipmissing(site["SALT"])))
    hideydecorations!(ax4, grid=false)

    ax5 = Axis(fig[1, 5], xlabel="ρ′ (kg/m³)", ylabel="z (m)")
    lines!(ax5, ρ′, z, linewidth=3)
    xlims!(ax5, extrema(skipmissing(site["RHOAnoma"])))
    hideydecorations!(ax5, grid=false)

    linkyaxes!(ax1, ax2, ax3, ax4, ax5)

    Label(fig[0, :], title, fontsize=30)

    mkpath(dirname(filepath))
    prog = Progress(length(frames), desc="Animating $filepath...", showspeed=true)
    record(fig, filepath, frames; framerate) do n
        frame[] = n
        next!(prog)
    end

    return nothing
end

# site = gather_site_data()
# plot_site_forcing(site)
# animate_site_profiles(site)
