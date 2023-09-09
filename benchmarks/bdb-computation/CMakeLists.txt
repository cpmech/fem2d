set(BENCHES
    bdb_classical
    bdb_expanded_full
    bdb_expanded
)

foreach(var ${BENCHES})
    add_executable(${var} "${var}.cpp")
    target_compile_definitions(${var} PUBLIC USE_MKL)
    target_link_libraries(${var} PUBLIC MKL::MKL ${LACLIB_LIBS} fem2d)
endforeach(var)
