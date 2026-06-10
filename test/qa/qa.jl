using QuasiMonteCarlo, Aqua, Test
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(QuasiMonteCarlo)
    Aqua.test_ambiguities(QuasiMonteCarlo, recursive = false)
    # check_extras disabled: `Pkg` is in [extras]/[targets].test without a [compat] bound.
    # Tracked in https://github.com/SciML/QuasiMonteCarlo.jl/issues/151
    Aqua.test_deps_compat(QuasiMonteCarlo, check_extras = false)
    @test_broken false  # Aqua deps_compat: `Pkg` extras dep lacks [compat] — tracked in https://github.com/SciML/QuasiMonteCarlo.jl/issues/151
    Aqua.test_piracies(QuasiMonteCarlo)
    Aqua.test_project_extras(QuasiMonteCarlo)
    Aqua.test_stale_deps(QuasiMonteCarlo)
    Aqua.test_unbound_args(QuasiMonteCarlo)
    Aqua.test_undefined_exports(QuasiMonteCarlo)
end
