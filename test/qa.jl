using QuasiMonteCarlo, Aqua
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(QuasiMonteCarlo)
    Aqua.test_ambiguities(QuasiMonteCarlo, recursive = false)
    Aqua.test_deps_compat(QuasiMonteCarlo)
    Aqua.test_piracies(QuasiMonteCarlo)
    Aqua.test_project_extras(QuasiMonteCarlo)
    Aqua.test_stale_deps(QuasiMonteCarlo)
    Aqua.test_unbound_args(QuasiMonteCarlo)
    Aqua.test_undefined_exports(QuasiMonteCarlo)
end
