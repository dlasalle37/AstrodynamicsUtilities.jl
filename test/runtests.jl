using AstrodynamicsUtilities
using Test

@testset "AstrodynamicsUtilities.jl" begin
    # Write your tests here.
end

@testset "State Representations" begin
    nothing_cart = CartesianState()
    @test isnothing(nothing_cart.vec)

    nothing_kep = KeplerianState()
    @test isnothing(nothing_kep.vec)

    # LEO test
    mu = 398600.4418
    ct = CartesianState([6678.0; 0.0; 0.0; 0.0; 7.07*sqrt(2)/2; 7.07*sqrt(2)/2])
    kp = transform(ct, KeplerianState(), mu)
    ctcheck = transform(kp, CartesianState(), mu)
    kpcheck = transform(ctcheck, KeplerianState(), mu)
    @test isapprox(ct.vec, ctcheck.vec; atol=1e-9)
    @test isapprox(kp.vec, kpcheck.vec; atol=1e-9)

    # GEO test
    mu = 398600.4418
    ct = CartesianState([42164; 0.0; 0.0; 0.0; 3.04*0.995037; 3.04*0.0995037])
    kp = transform(ct, KeplerianState(), mu)
    ctcheck = transform(kp, CartesianState(), mu)
    kpcheck = transform(ctcheck, KeplerianState(), mu)
    @test isapprox(ct.vec, ctcheck.vec; atol=1e-9)
    @test isapprox(kp.vec, kpcheck.vec; atol=1e-9)
    
end
