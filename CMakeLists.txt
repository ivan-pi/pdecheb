
cmake_minimum_required(VERSION 3.14)

project(pdecheb
    VERSION 1.0.0
    LANGUAGES Fortran
    DESCRIPTION "PDECHEB Package")

add_library(linpack linpack/linpack.f)

add_library(dassl
    dassl/ddaini.f
    dassl/ddajac.f
    dassl/ddanrm.f
    dassl/ddaslv.f
    dassl/ddassl.f
    dassl/ddastp.f
    dassl/ddatrp.f
    dassl/ddawts.f
    dassl/s88fmt.f
    dassl/scherr.f
    dassl/xerrwv.f)
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
    pdecheb/pdechb.f)

target_compile_options(pdecheb PRIVATE "$<$<Fortran_COMPILER_ID:GNU>:-std=legacy>")
