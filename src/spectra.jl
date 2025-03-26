using DataFrames
using FFTW
using FINUFFT

import DSP


"""
    nuwelch(x, timestamps, segment_length, Nf, overlap)

Perform non-uniform FFT with the Welch method.

https://en.wikipedia.org/wiki/Welch%27s_method
"""
function nuwelch(
        x::Array{<:Real},
        timestamps::Array{<:Dates.AbstractTime},
        segment_length::Day,
        Nf::Int;
        overlap::Real=0.5,
        )
    nuwelch(x .+ 0im, timestamps, segment_length, Nf; overlap)
end
export nuwelch


function nuwelch(
        x::Array{<:Complex},
        timestamps::Array{<:Dates.AbstractTime},
        segment_length::Day,
        Nf::Int;
        overlap::Real=0.5,
        )
    t = datetime2julian.(timestamps)
    t .-= t[1]
    n_segments = ceil(Int, (t[end] - t[1]) / segment_length.value / overlap) - 1
    coeffs = Array{Complex}(undef, n_segments, Nf)
    powers = Array{Real}(undef, n_segments, Nf)
    i = 1
    ti1 = 1
    ti2 = 1
    while ti2 < length(t)
        ti2 = min(searchsortedfirst(timestamps, timestamps[ti1] + segment_length), length(t))
        if ti2 - ti1 > floor(Int, 0.8 * segment_length.value * 8)
            tsegm = t[ti1: ti2]
            tsegm = tsegm .- minimum(tsegm)
            tsegm = tsegm * 2π / segment_length.value .- π
            window = DSP.hanning(ti2 - ti1 + 1)
            coeffs[i, :] = nufft1d1(tsegm, window .* x[ti1: ti2], -1, 1e-12, Nf, modeord=0)
            i += 1
        end
        ti1 = ti2 - floor(Int, (ti2 - ti1) * overlap)  # TODO: indices or duration?
    end
    coeffs = coeffs[1:i-1, :]
    powers = abs2.(coeffs)

    return mean(powers, dims=1)
end


"""
    welch(x, timestamps, segment_length, Nf, overlap)

Perform FFT with the Welch method.

https://en.wikipedia.org/wiki/Welch%27s_method
"""
function welch(
        x::Array{<:Complex},
        timestamps::Array{<:Dates.AbstractTime},
        segment_length::Day,
        Nf::Int;
        overlap::Real=0.5,
        )
    t = datetime2julian.(timestamps)
    t .-= t[1]
    n_segments = ceil(Int, (t[end] - t[1]) / segment_length.value / overlap) - 1
    coeffs = Array{Complex}(undef, n_segments, Nf)
    powers = Array{Real}(undef, n_segments, Nf)
    i = 1
    ti1 = 1
    ti2 = 1
    while ti2 < length(t)
        ti2 = min(searchsortedfirst(timestamps, timestamps[ti1] + segment_length), length(t))
        if (ti2 - ti1) == Nf
            window = DSP.hanning(ti2 - ti1)
            coeffs[i, :] = FFTW.fftshift(FFTW.fft(window .* x[ti1: ti2-1]))
            i += 1
        end
        ti1 = ti2 - floor(Int, Nf * overlap)
    end
    coeffs = coeffs[1:i-1, :]
    powers = abs2.(coeffs)

    return mean(powers, dims=1)
end
export welch
