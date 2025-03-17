using DataFrames
"""
Various calculations for preprocessing.
"""


"""
    haversine(ϕ₁, ϕ₂, λ₁, λ₂; R=6378100)

Calculate Haversine distance.

Coordinates should be given in degrees.

https://en.wikipedia.org/wiki/Haversine_formula
"""
function distance(lat_1::T, lat_2::T, lon_1::T, lon_2::T;
                   radius::T=6378.1e3) where T <: Real
    term1 = sind((lat_2 - lat_1) / 2.0)^2
    term2 = sind((lon_2 - lon_1) / 2.0)^2 * cosd(lat_1) * cosd(lat_2)
    d = 2.0 * radius * asin(sqrt(term1 + term2))
    return d
end


function distance(row1::DataFrameRow, row2::DataFrameRow)
    d = distance(row1[:lat], row2[:lat], row1[:lon], row2[:lon])
    return d
end


function distance(df::AbstractDataFrame)
    d = [distance(df[i, :], df[i+1, :]) for i=1:size(df, 1)-1]
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
    df[!, :distance] = distance(df0)
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
