using CairoMakie
using DataFrames
using DataFramesMeta
using GeoMakie

import ArchGDAL
import GADM


"Coastline for Bay of Bothnia"
function bob_coast()
    fin0 = GADM.get("FIN"; depth=0) |> DataFrame;
    swe0 = GADM.get("SWE"; depth=0) |> DataFrame;
    finswe = ArchGDAL.union(fin0.geom[1], swe0.geom[1])
    bobbounds = ArchGDAL.createmultipolygon([[[(lon, lat) for lon in (19, 26) for lat in (60, 66)][[1,2,4,3,1]]]])
    return ArchGDAL.intersection(finswe, bobbounds)
end


function plot_trajectories(df::AbstractDataFrame; color=:blue)
    fig = Figure(size=(600, 500))
    ga = GeoAxis(fig[1, 1]; dest = "+proj=natearth", title = "")
    # xlims!(lonlims...)
    # ylims!(latlims...)
    poly!(ga, bob_coast(), color=:grey)
    plt = scatter!(ga, df.lon, df.lat, color=color, markersize=5, alpha=1)
    # Colorbar(fig[1, 2], plt, label="Drift speed [m/s]")
    # df2 = combine(groupby(df, [^(:JP), ^(:winter)]), first)
    # scatter!(ga, df2.lon, df2.lat, color=^(:tomato), markersize=10, alpha=1)

    return fig
end
export plot_trajectories


function plot_trajectories2(df::AbstractDataFrame)
    fig = Figure(size=(600, 500))
    ga = GeoAxis(fig[1, 1]; dest = "+proj=natearth", title = "")
    # xlims!(lonlims...)
    # ylims!(latlims...)
    poly!(ga, bob_coast, color=^(:grey))
    plt = scatter!(ga, :lon, :lat, color=:speed, markersize=5, alpha=1)
    Colorbar(fig[1, 2], plt, label="Drift speed [m/s]")
    df2 = combine(groupby(df, [^(:JP), ^(:winter)]), first)
    scatter!(ga, df2.lon, df2.lat, color=^(:tomato), markersize=10, alpha=1)

    # ga2 = GeoAxis(fig; dest="+proj=natearth", bbox=BBox(340, 470, 80, 210))
    ga2 = Axis(fig; bbox=BBox(340, 470, 80, 210))
    hidedecorations!(ga2)
    # xlims2 = (-12, 40)
    # ylims2 = (30, 75)
    xlims2 = (4, 33)
    ylims2 = (53, 72)
    xlims!(ga2, xlims2...)
    ylims!(ga2, ylims2...)
    limbox = ArchGDAL.createmultipolygon([[[(lon, lat) for lon in xlims2 for lat in ylims2][[1,2,4,3,1]]]])
    subplt = poly!(ga2, limbox, strokecolor=^(:black), color=^(:white), strokewidth=8, overdraw=true, transparency=false)
    # subplt2 = lines!(ga2, GeoMakie.coastlines(), color=^(:black))
    subplt2 = image!(ga2, -180..180, -90..90, rotr90(GeoMakie.earth()), interpolate=false)
    bobbox = ArchGDAL.createmultipolygon([[[(lon, lat) for lon in lonlims for lat in latlims][[1,2,4,3,1]]]])
    subplt3 = poly!(ga2, bobbox, strokecolor=^(:red), strokewidth=3, color=^(:transparent))
    translate!(subplt, 0, 0, 100)
    translate!(subplt2, 0, 0, 101)
    translate!(subplt3, 0, 0, 102)

    return fig
end
