set(TESTS
    z_test_truss2d
)

foreach(var ${TESTS})
    add_executable(${var} "${var}.cpp")
    add_test(NAME ${var} COMMAND ${var})
    target_compile_definitions(${var} PUBLIC USE_MKL)
    target_link_libraries(${var} PUBLIC MKL::MKL ${LACLIB_LIBS} fem2d)
endforeach(var)
