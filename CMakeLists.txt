
cmake_minimum_required(VERSION 3.14)

project(pdecheb
    VERSION 1.0.0
    LANGUAGES Fortran
    DESCRIPTION "PDECHEB Package")

add_library(linpack linpack/linpack.f)

file(DOWNLOAD https://netlib.org/ode/ddassl.f ddassl.f
    EXPECTED_MD5 0c88cc5a08ecc31162bf0e3fead9365d)

add_library(dassl ddassl.f)

target_link_libraries(dassl PRIVATE linpack)
target_compile_options(dassl PRIVATE "$<$<Fortran_COMPILER_ID:GNU>:-std=legacy>")

add_library(pdecheb
    pdecheb/chintr.f
    pdecheb/cres.f
    pdecheb/cset.f
    pdecheb/dres.f
    pdecheb/error.f
    pdecheb/inichb.f
    pdecheb/interc.f
    pdecheb/intrch.f
    pdecheb/pdechb.f
    pdecheb/scherr.f)

target_compile_options(pdecheb PRIVATE "$<$<Fortran_COMPILER_ID:GNU>:-std=legacy>")

add_executable(ex1 examples/example1.f examples/dgejac.f)
target_link_libraries(ex1 pdecheb dassl)
target_compile_options(ex1 PRIVATE "$<$<Fortran_COMPILER_ID:GNU>:-std=legacy>")

add_executable(ex2 examples/example2.f examples/dgejac.f)
target_link_libraries(ex2 pdecheb dassl)
target_compile_options(ex2 PRIVATE "$<$<Fortran_COMPILER_ID:GNU>:-std=legacy>")

add_executable(ex3 examples/example3.f examples/dgejac.f)
target_link_libraries(ex3 pdecheb dassl)
target_compile_options(ex3 PRIVATE "$<$<Fortran_COMPILER_ID:GNU>:-std=legacy>")


add_executable(ex4 examples/example4.f examples/dgejac.f)
target_link_libraries(ex4 pdecheb dassl)
target_compile_options(ex4 PRIVATE "$<$<Fortran_COMPILER_ID:GNU>:-std=legacy>")

add_executable(ex5 examples/example5.f examples/dgejac.f)
target_link_libraries(ex5 pdecheb dassl)
target_compile_options(ex5 PRIVATE "$<$<Fortran_COMPILER_ID:GNU>:-std=legacy>")
