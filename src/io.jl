"""
Various utils for reading in and preprocessing ice drifter buoy data.

Warning: Very ad hoc.
"""

using CSV
using DataFrames
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


function get_icechart_dataset(year::Int)
    dsfile = DATAPATH * "icecharts/icecharts$(year).nc"
    ds = NCDataset(dsfile, "r")
    return ds
end
