using Test
using TJLF

# Regression for the SAT3 ITG/TEM discriminator NaN.
#
# intensity_sat picks x = |W_e| / sum_ions(W_i) at kmax by interpolating a ratio array. The zero
# guard brackets kmax with a strict `<` search (the pair *ending* at kmax) while
# linear_interpolation searches with `<=` (the pair *starting* at kmax). get_zonal_mixing returns
# kmax exactly on a grid point in four of its branches, and there the two brackets differ by one.
# When that extra point is a stable ky with zero ion QL weight, the whole-array ratio is 0/0 = NaN;
# the NaN's coefficient cancels algebraically at a node but not in floating point, so x came out
# NaN. QLA_P/QLA_E are then NaN and multiply every ky and every species, so the particle and energy
# fluxes of all species went NaN while the stress channels (QLA_O, a constant) survived.
#
# The three decks below are ITER local-equilibrium points that reproduced it (0.83% of a 30k-record
# ITER set did). No Fortran baseline is needed: the assertion is that the fluxes are finite, and
# that the survivor pattern (finite stress, NaN particle/energy) cannot reappear.
@testset "SAT3 x-discriminator finiteness" begin
    dir = joinpath(@__DIR__, "sat3_nan_regression")
    decks = sort([joinpath(dir, f) for f in readdir(dir) if startswith(f, "input.tglf.")])
    @test !isempty(decks)
    for deck in decks
        inputTJLF = readInput(deck)
        @test inputTJLF.SAT_RULE == 3
        out = TJLF.run(inputTJLF)
        @testset "$(basename(deck))" begin
            @test all(isfinite, out.QL_flux_out)
            @test all(isfinite, out.flux_spectrum)
            for (name, f) in (("Qe", TJLF.Qe), ("Qi", TJLF.Qi), ("Ge", TJLF.Γe), ("Pi", TJLF.Πi))
                @test isfinite(f(out.QL_flux_out))
            end
        end
    end
end
