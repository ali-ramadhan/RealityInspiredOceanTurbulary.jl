using ProgressMeter
using CairoMakie

function animate_les_profiles(ds; resolution=(1920, 1080), framerate=15, filepath=joinpath(config[:directories][:figures_dir], "les_profiles.mp4"))
    grid = ds["U"].grid
    times = ds["U"].times

    frames = 1:length(times)
    frame = Observable(1)

    title = @lift "LES profiles: $(times[$frame])"

    z = znodes(Center, grid)[:]
    U = @lift interior(ds["U"])[1, 1, :, $frame]
    V = @lift interior(ds["V"])[1, 1, :, $frame]
    T = @lift interior(ds["T"])[1, 1, :, $frame]
    S = @lift interior(ds["S"])[1, 1, :, $frame]
    B = @lift interior(ds["B"])[1, 1, :, $frame]

    fig = Figure(; resolution)

    ax1 = Axis(fig[1, 1], xlabel="U (m/s)", ylabel="z (m)")
    lines!(ax1, U, z, linewidth=3)
    # lines!(ax1, U_geo, z, linewidth=3, linestyle=:dot)
    xlims!(ax1, extrema(ds["U"]))

    ax2 = Axis(fig[1, 2], xlabel="V (m/s)", ylabel="z (m)")
    lines!(ax2, V, z, linewidth=3)
    # lines!(ax2, V_geo, z, linewidth=3, linestyle=:dot)
    xlims!(ax2, extrema(ds["V"]))
    hideydecorations!(ax2, grid=false)

    ax3 = Axis(fig[1, 3], xlabel="T (°C)", ylabel="z (m)")
    lines!(ax3, T, z, linewidth=3)
    xlims!(ax3, extrema(ds["T"]))
    hideydecorations!(ax3, grid=false)

    ax4 = Axis(fig[1, 4], xlabel="S (psu)", ylabel="z (m)")
    lines!(ax4, S, z, linewidth=3)
    xlims!(ax4, extrema(ds["S"]))
    hideydecorations!(ax4, grid=false)

    ax5 = Axis(fig[1, 5], xlabel="B (m/s²)", ylabel="z (m)")
    lines!(ax5, B, z, linewidth=3)
    xlims!(ax5, extrema(ds["B"]))
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
