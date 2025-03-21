using CairoMakie
using DataFrames
using DataFramesMeta
using GeoMakie

import ArchGDAL
import GADM


const LATLIMS = (63, 66)
const LONLIMS = (21, 25.1)


"""Create polygon for the coastline of Bay of Bothnia."""
function bob_coast()
    fin0 = GADM.get("FIN"; depth=0) |> DataFrame;
    swe0 = GADM.get("SWE"; depth=0) |> DataFrame;
    finswe = ArchGDAL.union(fin0.geom[1], swe0.geom[1])
    bobbounds = ArchGDAL.createmultipolygon([[[(lon, lat)
        for lon in (19, 26) for lat in (60, 66)][[1,2,4,3,1]]]])
    return ArchGDAL.intersection(finswe, bobbounds)
end


"""Create a background figure with map of BoB."""
function plot_bob_coast(;figsize=(600, 500), xlims=LONLIMS, ylims=LATLIMS)
    fig = Figure(size=figsize)
    ga = GeoAxis(fig[1, 1]; dest = "+proj=natearth", title = "")
    xlims!(xlims...)
    ylims!(ylims...)
    poly!(ga, bob_coast(), color=:grey)
    return (fig, ga)
end
export plot_bob_coast


"""Create a trajectory plot."""
function plot_trajectories(df::AbstractDataFrame; color=:blue, return_all=false)
    fig, ga = plot_bob_coast()
    plt = scatter!(ga, df.lon, df.lat, color=color, markersize=5, alpha=1)

    if return_all
        return (fig, ga, plt)
    else
    return fig
end
end
export plot_trajectories


"""Add box with map of northern Europe; BoB highlighted."""
function plot_add_insetmap!(fig::Figure, bbox::Rect2)
    ga2 = Axis(fig; bbox=bbox)
    hidedecorations!(ga2)
    xlims2 = (4, 33)  # Northern Europe
    ylims2 = (53, 72)
    xlims!(ga2, xlims2...)
    ylims!(ga2, ylims2...)
    limbox = ArchGDAL.createmultipolygon([[[(lon, lat)
        for lon in xlims2 for lat in ylims2][[1,2,4,3,1]]]])
    subplt = poly!(ga2, limbox, strokecolor=:black, color=:white,
        strokewidth=8, overdraw=true, transparency=false)
    subplt2 = image!(ga2, -180..180, -90..90, rotr90(GeoMakie.earth()),
        interpolate=false)
    bobbox = ArchGDAL.createmultipolygon([[[(lon, lat)
        for lon in LONLIMS for lat in LATLIMS][[1,2,4,3,1]]]])
    subplt3 = poly!(ga2, bobbox, strokecolor=:red, strokewidth=3, color=:transparent)
    translate!(subplt, 0, 0, 100)
    translate!(subplt2, 0, 0, 101)
    translate!(subplt3, 0, 0, 102)
    return fig
end
export plot_add_insetmap!


    return fig
end
