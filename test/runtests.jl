using Dates
using IceDrifters
using Test


@testset "io.jl" begin
    @test IceDrifters.parse_degrees("10d 30") == 10.5
    @test IceDrifters.get_jpnum("JP15") == 15
end


@testset "utils.jl" begin
    lons = [1.0, 2.0, 3.0, 4.0, 5.0]
    lats = [10.0, 9.0, 8.0, 7.0]
    times = [DateTime(2000, 1, d, h, 0) for d=1:2 for h = 0: 6: 23]
    vs = [stack([reshape(collect(1:20.0), 5, 4) .* 10^i for i = 0:7], dims=3)]
    @test interp_val(vs, lons, lats, times, 3.0, 8.0,
                     DateTime(2000, 1, 1, 0, 0))[1] ≈ 13.0
    @test interp_val(vs, lons, reverse(lats), times, 3.0, 8.0,
                     DateTime(2000, 1, 1, 0, 0))[1] ≈ 8.0
    @test interp_val(vs, lons, lats, Date.(times), 3.0, 8.0,
                     Date(2000, 1, 1))[1] ≈ (13.0 + 130 + 1300 + 13000) / 4
    @test interp_val(vs, lons, lats, times, 3.0, 8.0,
                     DateTime(2000, 1, 1, 3, 0))[1] ≈ (13.0 + 130) / 2
    @test interp_val(vs, lons, lats, times, 3.0, 8.5,
                     DateTime(2000, 1, 1, 0, 0))[1] ≈ (13.0 + 8.0) / 2
    @test interp_val(vs, lons, lats, times, 2.2, 8.0,
                     DateTime(2000, 1, 1, 0, 0))[1] ≈ 12.0 * 0.8 + 13.0 * 0.2
    @test interp_val(vs, lons, reverse(lats), times, 3.0, 8.2,
                     DateTime(2000, 1, 1, 0, 0))[1] ≈ 8.0 * 0.8 + 13.0 * 0.2
    @test interp_val(vs, lons, lats, times, 1.5, 9.2,
                     DateTime(2000, 1, 1, 0, 0))[1] ≈ ((1 + 2) * 0.2 + (6 + 7) * 0.8) / 2

    @test closest_val(vs, lons, reverse(lats), times, 1.5, 9.2,
                     DateTime(2000, 1, 1, 0, 0))[1] ≈ 12.0
    @test closest_val(vs, lons, lats, times, 3.5, 9.2,
                     DateTime(2000, 1, 1, 0, 0))[1] ≈ 9.0
    @test closest_val(vs, lons, lats, times, 1.8, 9.2,
                     DateTime(2000, 1, 1, 10, 0))[1] ≈ 700.0

    @test winter(DateTime(2014, 2, 1, 10, 0)) == 2014
    @test winter(DateTime(2014, 9, 1, 10, 0)) == 2015
end


@testset "buoy_calculations.jl" begin
    @test IceDrifters.distance(-20., 20., 0., 180.) ≈ π * 6378100
    @test IceDrifters.bearing(0, 10, 0, 0) ≈ 0
    @test IceDrifters.bearing(10, 0, 0, 0) ≈ 180
    @test IceDrifters.bearing(0., 0.1, 0., 0.1) ≈ 45 atol=1e-3
    @test IceDrifters.bearing(0., -0.1, 0., -0.1) ≈ 225 atol=1e-3
    @test IceDrifters.uv2dir(1, 1) ≈ 45 atol=1e-3
    @test IceDrifters.uv2dir(-1, -1) ≈ 225 atol=1e-3
end
