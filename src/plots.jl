using CairoMakie
using DataFrames
using DataFramesMeta
using GeoMakie

import ArchGDAL
import GADM


const LONLIMS = (21, 25.1)
const LATLIMS = (63, 66)
export LONLIMS, LATLIMS


"""Add some noise to data (for plotting of categories)."""
jitter(data; scale=1.0) = data + (rand(length(data)) .- 0.5) * scale
export jitter


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


function plot_triangle(t::Triangle; labels=true)
    fig = Figure()
    ax = Axis(fig[1,1])
    lines!(ax, [t.p1.x, t.p2.x], [t.p1.y, t.p2.y])
    labels && text!(ax, t.p1.x, t.p1.y, text="1")
    lines!(ax, [t.p2.x, t.p3.x], [t.p2.y, t.p3.y])
    labels && text!(ax, t.p2.x, t.p2.y, text="2")
    lines!(ax, [t.p3.x, t.p1.x], [t.p3.y, t.p1.y])
    labels && text!(ax, t.p3.x, t.p3.y, text="3")
    return fig
end
export plot_triangle


function plot_triangles(df::AbstractDataFrame)
    fig, ga = plot_bob_coast()
    plot_triangles!(df, ga)
    return fig
end
export plot_triangles


function plot_triangles!(df::AbstractDataFrame, ga::GeoAxis)
    for points in gather_latlons(df)
        scatterlines!(ga, points, markercolor=:black, markersize=15)
    end
end
export plot_triangles!


function animate_triangles!(gdf::GroupedDataFrame; framerate::Int=5)
    fig, ga = plot_bob_coast()
    plots = []
    rec = Record(fig, gdf, framerate=framerate) do df
        while length(plots) > 0
            delete!(ga, pop!(plots))
        end
        latlons = gather_latlons(df)
        for points in latlons
            push!(plots, scatterlines!(ga, points))
        end
        ga.title = string(Date(df[1, :time]))
    end

    return rec
end
export animate_triangles!


function plot_deform_split(df::AbstractDataFrame, splitcol::Symbol, splitval::Real;
        color=to_colormap(:tab10)[1], kwargs...)
    df = @rsubset(df, :deform > 0)
    fig = Figure(size=(700, 500))
    splitcolstr = replace(uppercasefirst(string(splitcol)), "_"=>" ")
    i = 1
    sel = df[!, splitcol] .> splitval
    println("Split: $(sum(sel)) / $(nrow(df))")
    for sel1 in (.!sel, sel)
        df_ = df[sel1, :]
        scale = df_.scale / 1000
        ax = Axis(fig[1, i], xscale=log10, yscale=log10,
            xlabel="Length scale [km]", ylabel="Deformation [1/day]",
            xtickformat="{:d}", xticks=[1,10,100])
        (i == 1) ? (ax.title = "a) $splitcolstr ≤ $splitval") :
                   (ax.title = "b) $splitcolstr > $splitval")
        if color isa Symbol
            colori = getproperty(df_, color)
        else
            colori = color
        end
        scatter!(ax, scale, df_.deform; alpha=0.9, markersize=5, color=colori, kwargs...)
        pl = get_power_law(df_)
        lines!(ax, scale, predict_power_law.(scale, Ref(pl)), color=:red,
            label=L"D=%$(round(pl.α, digits=1))L^{%$(round(pl.β, digits=2))}")
        axislegend(ax, position=:rt, framevisible=false)
        xlims!(ax, 6, 140)
        ylims!(ax, 1e-3, 1e2)
        (i == 2) && hideydecorations!(ax, grid=false)
        i += 1
    end
    linkaxes!([a for a in fig.content if typeof(a) == Axis]...)
    return fig
end
export plot_deform_split
