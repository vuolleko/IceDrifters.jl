using DataFrames
using NCDatasets
using Statistics


"""
    interp_val(vs, lon, lat, t, ds)

Interpolate to find a value for one or more symbols in a NetCDF file.
"""
function interp_val(vs::Vector{Symbol}, lon::Real, lat::Real, t::TimeType, ds::NCDataset)
    lons = view(ds[:longitude], :)  # ascending
    lats = view(ds[:latitude], :)  # descending
    times = view(ds[:valid_time], :)
    arr = [view(ds[v], :, :, :) for v in vs]
    return interp_val(arr, lons, lats, times, lon, lat, t)
end
export interp_val


function interp_val(
        vs::Vector{<:AbstractArray{T, 3}},
        lons::AbstractArray{<:Real},
        lats::AbstractArray{<:Real},
        times::AbstractArray{<:TimeType},
        lon::Real,
        lat::Real,
        t::TimeType,
        ) where T <: Union{Missing, Real}
    i = searchsortedfirst(lons, lon, rev=(lons[1] > lons[2]))
    j = searchsortedfirst(lats, lat, rev=(lats[1] > lats[2]))
    ks = searchsorted(times, t)  # note: removed casting to same type
    (ks.start > ks.stop) && (ks = ks.stop: ks.start)  # case when exact match not found
    dlon1 = (lons[i] - lon) / (lons[i] - lons[i-1])
    dlon2 = (lon - lons[i-1]) / (lons[i] - lons[i-1])
    dlat1 = (lats[j-1] - lat) / (lats[j-1] - lats[j])
    dlat2 = (lat - lats[j]) / (lats[j-1] - lats[j])
    vals = Vector{Float64}(undef, length(vs))
    for (ii, v) in enumerate(vs)
        vals1 = dlat1 * (dlon1 * v[i-1, j, ks] + dlon2 * v[i, j, ks]) +
                dlat2 * (dlon1 * v[i-1, j-1, ks] + dlon2 * v[i, j-1, ks])
        vals[ii] = mean(vals1)
    end
    return vals
end


interp_val(v::Symbol, lon::Real, lat::Real, t::TimeType, ds::NCDataset) =
    interp_val([v], lon, lat, t, ds)[1]


"""
    closest_val(vs, lon, lat, t, ds)

Find the closest value for symbols in a NetCDF file.
"""
function closest_val(vs::Vector{Symbol}, lon::Real, lat::Real, t::TimeType, ds::NCDataset)
    lats = ds[:lat][:]
    lons = ds[:lon][:]
    times = ds[:TIME][:]
    arr = [view(ds[v], :, :, :) for v in vs]
    return closest_val(arr, lons, lats, times, lon, lat, t)
end
export closest_val


function closest_val(
        vs::Vector{<:AbstractArray{T, 3}},
        lons::Vector{<:Real},
        lats::Vector{<:Real},
        times::Vector{<:TimeType},
        lon::Real,
        lat::Real,
        t::TimeType,
        ) where T <: Union{Missing, Real}
    i = searchsortedfirst(lons, lon, rev=(lons[1] > lons[2]))
    j = searchsortedfirst(lats, lat, rev=(lats[1] > lats[2]))
    abs(lons[i] - lon) > abs(lon - lons[i-1]) && (i -= 1)
    abs(lats[j] - lat) > abs(lat - lats[j-1]) && (j -= 1)
    k = searchsortedfirst(typeof(t).(times), t)
    vals = Vector{Union{Float64, Missing}}(undef, length(vs))
    for (ii, v) in enumerate(vs)
        vals[ii] = v[i, j, k]
    end
    return vals
end


closest_val(cols::Tuple, row::DataFrameRow, ds::NCDataset) =
    closest_val([cols...], row[:lon], row[:lat], Date(row[:timestamp]), ds)


winter(dt::DateTime) = year(dt) + (month(dt) > 8)  # winter 2011-2012 marked as 2012 etc.
export winter
