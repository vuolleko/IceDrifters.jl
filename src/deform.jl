"""
Deformation analysis; calculations on imaginary triangles formed by buoys.
"""

using Combinatorics
using DataFrames
using DataFramesMeta
using LinearAlgebra
using Statistics


struct Vertex
    x::Real
    y::Real
end
Vertex(r::DataFrameRow) = Vertex(r[:x], r[:y])


import Base:-
-(v1::Vertex, v2::Vertex) = [v1.x - v2.x, v1.y - v2.y]


struct Triangle
    p1::Vertex
    p2::Vertex
    p3::Vertex
end
Triangle(x1::T, y1::T, x2::T, y2::T, x3::T, y3::T) where
    T<:Real = Triangle(Vertex(x1, y1), Vertex(x2, y2), Vertex(x3, y3))
Triangle(df::AbstractDataFrame) = Triangle([Vertex(row) for row in eachrow(df)]...)


area(t::Triangle) = 0.5 * det([1 1 1; t.p1.x t.p2.x t.p3.x; t.p1.y t.p2.y t.p3.y])


isposorient(t::Triangle) = area(t) > 0


import LinearAlgebra: normalize
normalize(v::Vertex) = [v.x, v.y] ./ sqrt(v.x^2 + v.y^2)


function angles(t::Triangle)
    v1 = normalize(t.p2 - t.p1)
    v2 = normalize(t.p3 - t.p2)
    v3 = normalize(t.p1 - t.p3)
    return acosd.(clamp.([-v1' * v2, -v2' * v3, -v3' * v1], -1, 1))
end


istoosharp(t::Triangle, α::Number) = minimum(angles(t)) < α


vertices(t::Triangle) = [t.p1.x, t.p1.y, t.p2.x, t.p2.y, t.p3.x, t.p3.y]


"""
    get_triangles(df; df_static, max_static, keep_smallest_static)

Form all eligible triangles and calculate their differential kinematic properties.

Triangles are formed from buoys with matching timestamp. Result sorted by time.

Arguments
=========
- df: DataFrame with buoy data.
- df_static: DataFrame of fixed reference buoys.
- max_static: Maximum number of fixed reference buoys in a triangle.
- keep_smallest_static: keep only the smallest triangle, when static components.
"""
function get_triangles(
        df::AbstractDataFrame;
        df_static::Union{AbstractDataFrame, Nothing}=nothing,
        max_static::Int=1,
        keep_smallest_static::Bool=true,
        )::DataFrame

    df2 = transform(df, [:lon, :lat] => ByRow(coordtrans()) => [:x, :y])
    gdf = groupby(df2, :timestamp)
    triangles = get_triangles(gdf, df_static=df_static, max_static=max_static,
        keep_smallest_static=keep_smallest_static)
    return sort(triangles, :time)
end
export get_triangles


"""
    get_triangles(gdf; df_static, max_static, keep_smallest_static)

Form all eligible triangles and calculate their differential kinematic properties.

Arguments
=========
- gdf: grouped DataFrame e.g. by timestamp.
- df_static: DataFrame of fixed reference buoys.
- max_static: Maximum number of fixed reference buoys in a triangle.
- keep_smallest_static: keep only the smallest triangle, when static components.
"""
function get_triangles(
        gdf::GroupedDataFrame;
        df_static::Union{AbstractDataFrame, Nothing}=nothing,
        max_static::Int=1,
        keep_smallest_static::Bool=true,
        )::DataFrame

    triangles = DataFrame(
        time=DateTime[], x1=Real[], y1=Real[], x2=Real[], y2=Real[],
        x3=Real[], y3=Real[], area=Real[], u1=Real[], u2=Real[], u3=Real[],
        v1=Real[], v2=Real[], v3=Real[], JP=Vector{<:Int}[],
        )
    any_static = !isnothing(df_static)
    n_static = any_static ? nrow(df_static) : 0

    # TODO: split into smaller functions
    for df in gdf
        # add static points
        any_static && (df = vcat(df, df_static, cols=:union))
        n = size(df, 1)
        if size(df, 1) < 3
            continue
        end

        n_combs = 0
        combs = combinations(1:size(df, 1), 3)
        for comb in combs
            # reject too static combinations
            any_static && (mapreduce(i->i in (n-n_static+1):n, +, comb) > max_static) && continue

            t = Triangle(df[comb, :])
            if ~isposorient(t)  # switch orientation
                comb = [comb[1], comb[3], comb[2]]
                t = Triangle(df[comb, :])
            end
            istoosharp(t, 15) && continue  # reject if any angle < 15 degrees
            u1, u2, u3 = df[comb, :speed] .* sind.(df[comb, :bearing])  # bearing 0 is north
            v1, v2, v3 = df[comb, :speed] .* cosd.(df[comb, :bearing])
            push!(triangles, [df[1, :timestamp], vertices(t)..., area(t),
                u1, u2, u3, v1, v2, v3, df[comb, :JP]])
            n_combs += 1
        end

        # only keep the smallest triangle with static
        if any_static & keep_smallest_static & (n_combs > 1)
            inds = findall(jps -> 0 in jps, triangles[end-n_combs+1: end, :JP])
            inds .+= size(triangles)[1] - n_combs
            isempty(inds) && continue
            df_ = triangles[inds, :]
            dfmin_ = combine(
                groupby(transform(df_, :JP => ByRow(Set) => :JPset), :JPset),
                :area => minimum)
            delinds = filter(i -> triangles[i, :area] ∉ dfmin_.area_minimum, inds)
            deleteat!(triangles, delinds)
        end
    end
    add_diffkinprops!(triangles)

    return triangles
end


function add_diffkinprops!(df::AbstractDataFrame)
    # calculate strain rates [1/s]
    @rtransform!(df, :dudx = 1/(2*:area) * ((:u1+:u3) * (:y1-:y3)
        + (:u1+:u2) * (:y1-:y2) + (:u2+:u3) * (:y2-:y3)))
    @rtransform!(df, :dudy = -1/(2*:area) * ((:u1+:u3) * (:x1-:x3)
        + (:u1+:u2) * (:x1-:x2) + (:u2+:u3) * (:x2-:x3)))
    @rtransform!(df, :dvdx = 1/(2*:area) * ((:v1+:v3) * (:y1-:y3)
        + (:v1+:v2) * (:y1-:y2) + (:v2+:v3) * (:y2-:y3)))
    @rtransform!(df, :dvdy = -1/(2*:area) * ((:v1+:v3) * (:x1-:x3)
        + (:v1+:v2) * (:x1-:x2) + (:v2+:v3) * (:x2-:x3)))

    # calculate divergence [1/day]
    @rtransform!(df, :div = (:dudx + :dvdy) * 24 * 3600)

    # calculate shear rate [1/day]
    @rtransform!(df, :shr = sqrt((:dudx - :dvdy)^2 + (:dudy + :dvdx)^2) * 24 * 3600)

    # calculate total deformation [1/day]
    @rtransform!(df, :deform = sqrt(:div^2 + :shr^2))
end


"""
    get_triangle_rows(trow, df)

Get subset of buoy data in df that relates to a triangle.
"""
get_triangle_rows(trow::DataFrameRow, df::AbstractDataFrame) = subset(df,
    :timestamp => t -> t .== trow.time, :JP => ByRow(j -> j in trow.JP))
export get_triangle_rows


"""
    add_triangle_mean_cols(tri_df, df)

Add several variables as mean for buoy locations in a triangle.
"""
function add_triangle_mean_cols!(tri_df::AbstractDataFrame, df::AbstractDataFrame)
    n = nrow(tri_df)
    newcols = [:dist2shore, :t_air, :wind_speed, :speed, :sithic, :siconc, :lon, :lat]
    for col in newcols
        tri_df[!, col] = Vector{Float64}(undef, n)
    end

    Threads.@threads for i=1:n
        df1 = get_triangle_rows(tri_df[i, :], df)
        for col in newcols
            tri_df[i, col] = mean(df1[!, col])
        end
    end
    transform!(tri_df,
        [:speed, :wind_speed] => ByRow((vi, vw) -> vi / vw) => :wind_factor)
    transform!(tri_df, :area => ByRow(sqrt) => :scale)
    return
end
export add_triangle_mean_cols!


function gather_latlons(df::AbstractDataFrame)
    latlons = Vector{Tuple{Real, Real}}[]
    for row in eachrow(df)
        points = map(inv(coordtrans()), zip(row[[:x1, :x2, :x3]], row[[:y1, :y2, :y3]]))
        push!(points, points[1])  # close the triangle
        push!(latlons, points)
    end
    return latlons
end


struct PowerLaw
    α::Real
    β::Real
end


"""
    get_power_law(df)

Fit a power law for scale vs. deformation rate.
"""
function get_power_law(df::AbstractDataFrame)
    @assert :deform in propertynames(df)
    df = subset(df, :deform => d -> d .> 0)
    scale = sqrt.(df.area) / 1000
    log_α, β = ([ones(length(scale)) log.(scale)] \ log.(df.deform))'
    return PowerLaw(exp(log_α), β)
end
export get_power_law


"""
    predict_power_law(scale, pl)

Predict deformation rate for given scale and power law.
"""
predict_power_law(scale::Real, pl::PowerLaw) = pl.α * scale ^ pl.β
export predict_power_law
