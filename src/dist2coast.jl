using GeoJSON
using SphericalGeometry


function create_polygon(geojsonfile::AbstractString, lon_inside::Real, lat_inside::Real)
    geojson = GeoJSON.read(geojsonfile)
    pts = SphericalGeometry.Point{Float64}[]
    for p in geojson
        push!(pts, SphericalGeometry.Point(Float64.(p)...))
    end

    inside = SphericalGeometry.Point(lon_inside, lat_inside)
    polygon = Polygon(inside, pts)
    return polygon
end


function dist2coast(lon::Real, lat::Real, polygon::Polygon)
    d = angular_distance(SphericalGeometry.Point(lon, lat), polygon)
    d = SphericalGeometry.distance(d)
    return d
end


dist2coast(lon::Real, lat::Real) = dist2coast([lon], [lat])[1]


const BOB_GEOJSON = "geojson/bob10.geojson"
const HAILUOTO_GEOJSON = "geojson/hailuoto.geojson"
const BOB_INSIDE = (23., 65.)
const HAILUOTO_INSIDE = (24.71, 65.03)


"""
    dist2coast(lons, lats)

Calculate distance to coast in Bay of Bothnia.
"""
function dist2coast(lons::AbstractArray{<:Real}, lats::AbstractArray{<:Real})
    bob = create_polygon(BOB_GEOJSON, BOB_INSIDE...)
    hailuoto = create_polygon(HAILUOTO_GEOJSON, HAILUOTO_INSIDE...)
    n = length(lons)
    dists = Vector{Float64}(undef, n)

    Threads.@threads for i = 1:n
        d_bob = dist2coast(lons[i], lats[i], bob)
        d_hailuoto = dist2coast(lons[i], lats[i], hailuoto)
        dists[i] = min(d_bob, d_hailuoto)
    end
    return dists
end
export dist2coast


isenclosed(lon::Real, lat::Real) = isenclosed([lon], [lat])[1]


"""
    isenclosed(lons, lats)

Check whether coordinates lie within the Bay of Bothnia sea area.
"""
function isenclosed(lons::AbstractArray{<:Real}, lats::AbstractArray{<:Real})
    n = length(lons)
    bob = create_polygon(BOB_GEOJSON, BOB_INSIDE...)
    arein = Vector{Bool}(undef, n)
    Threads.@threads for i = 1:n
        arein[i] = isinside(SphericalGeometry.Point(lons[i], lats[i]), bob)
    end
    return arein
end
export isenclosed
