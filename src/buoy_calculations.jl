"""
Various calculations for preprocessing.
"""


using DataFrames
using Dates


"""
    haversine(ϕ₁, ϕ₂, λ₁, λ₂; R=6378100)

Calculate Haversine distance.

Coordinates should be given in degrees.

https://en.wikipedia.org/wiki/Haversine_formula
"""
function get_distance(lat_1::T, lat_2::T, lon_1::T, lon_2::T;
                   radius::Real=6378.1e3) where T <: Real
    term1 = sind((lat_2 - lat_1) / 2.0)^2
    term2 = sind((lon_2 - lon_1) / 2.0)^2 * cosd(lat_1) * cosd(lat_2)
    d = 2.0 * radius * asin(sqrt(term1 + term2))
    return d
end


function get_distance(row1::DataFrameRow, row2::DataFrameRow)
    d = get_distance(row1[:lat], row2[:lat], row1[:lon], row2[:lon])
    return d
end


function get_distance(df::AbstractDataFrame)
    d = [get_distance(df[i, :], df[i+1, :]) for i=1:size(df, 1)-1]
    return d
end


"""Calculate bearing (north 0 degrees, east 90 degrees)."""
function bearing(lat_1::T, lat_2::T, lon_1::T, lon_2::T) where T <: Real
    Δlon = lon_2 - lon_1
    y = cosd(lat_2) * sind(Δlon)
    x = cosd(lat_1) * sind(lat_2) - sind(lat_1) * cosd(lat_2) * cosd(Δlon)
    b = atand(y, x)
    b = b < 0 ? b + 360 : b  # adjust [-180, 0[ -> [180, 360[
    return b
end


function bearing(df::AbstractDataFrame)
    return [bearing(df[i, :lat], df[i+1, :lat], df[i, :lon], df[i+1, :lon])
            for i = 1: nrow(df)-1]
end


"Calculate direction from velocity."
uv2dir(u, v) = mod(360 - (passmissing(atand)(v, u) - 90), 360)


"""
    distspeeds1(df)

Calculate travelled distance, speed and bearing for single buoy.
"""
function distspeeds1(df0::AbstractDataFrame; remove_gaps=false)
    @assert length(unique(df0.JP)) == 1
    @assert issorted(df0.timestamp)
    df = copy(df0[2:end, :])
    df[!, :distance] = get_distance(df0)
    dt = Second.(diff(df0.timestamp))
    df[!, :speed] = df[!, :distance] ./ Dates.value.(dt)  # [m/s]
    df[!, :bearing] = bearing(df0)

    if remove_gaps
        bad = dt .> Week(1)
        df = df[.!bad, :]
    end
    return df
end


"""
    distspeeds(df)

Calculate travelled distance, speed and bearing.
"""
function distspeeds(df0::AbstractDataFrame)
    df = transform(df0, :timestamp => ByRow(winter) => :winter)
    df2 = DataFrame()
    for df_ in groupby(df, [:winter, :JP])
        df1 = distspeeds1(df_)
        df2 = vcat(df2, df1)
    end
    return df2
end
export distspeeds


function average_3h(df0::AbstractDataFrame)
    df_tempo = transform(df0, :timestamp => (t -> hour.(t) .÷ 3) => :t3h);  # div
    meancols = [:sithic, :siconc, :t_air, :u_wind, :v_wind, :dist2shore, :sst]
    df_tempo = combine(
        groupby(df_tempo, [:doy, :t3h]),
        Not(meancols) .=> first .=> Not(meancols),
        meancols .=> mean .=> meancols,
        )
    return select(df_tempo, Not(:t3h))
end


function find_suspicious(df::AbstractDataFrame)
    @assert length(unique(df.JP)) == 1
    inds0 = findall(df.speed .== 0.0)
    isempty(inds0) && return Int[]
    (inds0[1] == 1) && (inds0 = inds0[2:end])
    (inds0[end] == nrow(df)) && (inds0 = inds0[1:end-1])
    inds_del = Int[]
    j = 1
    while j < length(inds0)
        prev = inds0[j] - 1
        while (j < length(inds0)) && (inds0[j+1] == inds0[j] + 1)
            j += 1
        end
        next = inds0[j] + 1
        (next == nrow(df)) && break
        if (min(df.speed[prev], df.speed[next+1]) > 0.02) && (df.speed[next] > (max(df.speed[prev], df.speed[next+1]) * 1.5))
            push!(inds_del, collect(prev+1: next-1)...)
        end
        j += 1
    end
    return inds_del
end


"""
    get_dts1(df)

Find distribution of timesteps between observations for one buoy.
"""
function get_dts1(df1::AbstractDataFrame)
    @assert length(unique(df1.JP)) == 1
    @assert length(unique(df1.winter)) == 1

    moves = cumsum(df1.distance .> 0)
    dts = Hour[]
    for m in unique(moves)
        ind1 = findfirst(moves .== m)
        ind2 = findlast(moves .== m) + 1
        (ind2 > nrow(df1)) && continue
        dt = df1[ind2, :timestamp] - df1[ind1, :timestamp]  # [ms]
        push!(dts, Hour(dt.value ÷ 3600000))
    end
    return sort(Dict(dt1 => sum(dt1 .== dts) / length(dts) for dt1 in unique(dts)))
end


"""
    get_dts(df)

Find distributions of timesteps between observed movement.

Some buoys update their position (GPS fix) less often than they report it.
Position for a given timestamp can be wrong. True update frequency can be
inferred from the distribution of timestemps between observed movement.
"""
function get_dts(df::AbstractDataFrame)
    dt_df = combine(
        groupby(df, [:winter, :JP]),
        get_dts1)
    transform!(dt_df, :x1 => ByRow(argmax) => :most_common)
    return dt_df
end
export get_dts


"""
    lonlat2gid(lon, lat)

 Assign each point into a 0.5x0.5 degree grid using an integer id.
 """
lonlat2gid(lon::Real, lat::Real) = floor.(Int, lat * 2) * 100 + floor.(Int, lon * 2)


"""
    gid2lonlat(x)

Convert the 0.5x0.5 degree cell id back to (lon, lat).
"""
function gid2lonlat(x::Int)
    xs = string(x)
    return (parse(Float64, xs[4:5]) / 2 + 0.25, parse(Float64, xs[1:3]) / 2 + 0.25)
end


struct Strides
    x::Vector{<:Real}  # momentary positions
    y::Vector{<:Real}
    sithic::Vector{<:Real}

    distances::Vector{<:Real}  # cumulative for each stride
    durations::Vector{Hour}
    doy::Vector{Int}  # final for each stride

    durations_all::Ref{Hour}  # for all movements
    strides::Dict{Int, Int}  # start ind, length

    Strides() = new(Real[], Real[], Real[], Real[], Hour[], Int[], Hour(0), Dict{Int, Int}())
end


"""
    get_strides(df; angle_max, speed_min, duration_min)

Perform persistence analysis on data.
"""
function get_strides(df::AbstractDataFrame;
                     angle_max::Real=30,
                     speed_min::Real=0.0,
                     duration_min::Hour=Hour(3),
                     )::Strides
    s = Strides()
    for jdf in groupby(df, [:JP, :winter])
        x1 = 0.0
        y1 = 0.0
        stride_ind1 = 1
        distance = 0.0
        duration = Hour(0)
        prev = first(jdf)
        for (ind, row) in enumerate(eachrow(jdf))
            Δθ = acosd(clamp(
                    (normalize([row.u row.v]) * normalize([prev.u, prev.v]))[1],
                    -1, 1))
            Δt = row.timestamp - prev.timestamp

            if (Δt > Day(1)) | (Δθ > angle_max) | (ind == nrow(jdf)) | (row.speed < speed_min)
                x1 = 0.0
                y1 = 0.0
                push!(s.x, NaN, 0)  # discontinuously back to 0
                push!(s.y, NaN, 0)
                push!(s.sithic, 1, 1)
                stridelen = ind - 1 - stride_ind1
                (stridelen > 8) && (s.strides[stride_ind1] = stridelen)
                stride_ind1 = ind - 1
                if duration >= duration_min  # record at stride end
                    push!(s.durations, duration)
                    push!(s.distances, distance)
                    push!(s.doy, row.doy)
                end
                duration = Hour(0)
                distance = 0.0

            else  # accumulate only if no reset
                duration += Hour(Δt - mod(Δt, 3600e3))
                distance += row.distance
                x1 += row.distance * sind(row.bearing)
                y1 += row.distance * cosd(row.bearing)
            end
            s.durations_all[] += Hour(Δt - mod(Δt, 3600e3))
            push!(s.x, x1)
            push!(s.y, y1)
            push!(s.sithic, row.sithic)
            prev = row
        end
    end
    return s
end
export get_strides
