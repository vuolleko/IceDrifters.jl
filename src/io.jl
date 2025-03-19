"""
Various utils for reading in and preprocessing ice drifter buoy data.

Warning: Very ad hoc.
"""

using CSV
using DataFrames
using DataFramesMeta
using Dates
using NCDatasets
using ZipArchives


function get_jp(num::Int)
    get_jp("JP" * string(num) * ".csv")
end


"""
    get_jp(filename)

Load buoy positions from raw CSV.
"""
function get_jp(filename::AbstractString)
    df0 = DataFrame(CSV.File(filename, dateformat="yyyy-mm-dd HH:MM:SS"))
    rename!(df0, Dict([s => strip(s) for s in names(df0)]))
    if "GPSDELAY" in names(df0)
        subset!(df0, :GPSDELAY => ByRow(d -> d < 6))  # allow 5 min delay
    elseif "TSLF" in names(df0)
        subset!(df0, :TSLF => ByRow(d -> d < 6))  # allow 5 min delay
    end
    df = select(df0,
        [:YEAR, :MONTH, :DAY, :HOUR, :MIN] =>
            ByRow((y, m, d, h, mi) -> DateTime(y, m, d, h, mi)) => :timestamp,
        :LONGITUDE => ByRow(parse_degrees) => :lon,
        :LATITUDE => ByRow(parse_degrees) => :lat,
        )
    df.sst = "SST" in names(df0) ? df0.SST : repeat([missing], nrow(df))
    df[!, :doy] = dayofyear.(df[!, :timestamp]);
    if "Asset Name" in names(df0)
        df.JP = map(get_jpnum, df0[!, "Asset Name"])
    else
        jp = split(filename, ['/', '.', '_'])[end-1]
        if startswith(jp, "JP")
            df.JP .= get_jpnum(jp)
        end
    end
    return df
end
export get_jp


get_jpnum(s::AbstractString) = parse(Int, s[3:end])


function parse_degrees(coord::AbstractString)
    d, m = split(coord)
    whole = parse(Float64, d[1:end-1])
    part = parse(Float64, m)
    return whole + part / 60.0
end


parse_degrees(coord::Real) = coord  # some CSVs have this in correct form


function get_zipped_icechart_dataset(year::Int)
    zipfile = DATAPATH * "icecharts/icecharts$(year).zip"
    arch = ZipReader(read(open(zipfile)))
    datavec = zip_readentry(arch, "icecharts$(year).nc")
    ds = NCDataset("dummy", "r", memory=datavec)
    return ds
end


"""
    get_era5_winds(df, ds)

Get interpolated winds and air temperature from a NetCDF file with ERA5 data.
"""
function get_era5_winds(df0::AbstractDataFrame, ds::NCDataset)
    df = copy(df0)
    t = Vector{Union{Float64, Missing}}(undef, nrow(df))
    u = Vector{Union{Float64, Missing}}(undef, nrow(df))
    v = Vector{Union{Float64, Missing}}(undef, nrow(df))
    for i = 1:nrow(df)
        (i % 5000 == 0) && println("$i/$(nrow(df)) done")
        row = view(df, i, :)
        t[i], u[i], v[i] = interp_val([:t, :u, :v], row[:lon], row[:lat], row[:timestamp], ds)
    end
    df[!, :t_air] = t .- 273.15;
    df[!, :u_wind] = u;
    df[!, :v_wind] = v;
    @rtransform!(df, :wind_speed = sqrt(:u_wind^2 + :v_wind^2))
    @rtransform!(df, :wind_direction = uv2dir(:u_wind, :v_wind))
    @rtransform!(df, :wind_factor = :speed / :wind_speed)
    return df
end
export get_era5_winds


const ICECHART_COLS = (;
    :sitype => :kFmiIceType,
    :siconc => :kFmiIceConcentration,
    :sithic => :kFmiIceThickness,
    :sithic_min => :kFmiIceMinThickness,
    :sithic_max => :kFmiIceMaxThickness,
    :siridging => :kFmiIceDegreeOfRidging,
    )


"""
    get_fmi_ice(df)

Get closest values for sea ice properties in ice charts.
"""
function get_fmi_ice(df0::AbstractDataFrame, filename_f::Function)
    df = transform(df0, :timestamp => ByRow(year) => :year)
    for k in keys(ICECHART_COLS)
        df[!, k] = Vector{Union{Float64, Missing}}(undef, nrow(df))
    end

    for d in groupby(df, :year)
        y = first(d)[:year]
        println("Processing icecharts for $y...")
        ds = NCDataset(filename_f(y), "r")
        for row in eachrow(d)
            row[[keys(ICECHART_COLS)...]] .= closest_val(values(ICECHART_COLS), row, ds)
        end
        close(ds)
    end

    return select(df, Not(:year))
end
export get_fmi_ice
