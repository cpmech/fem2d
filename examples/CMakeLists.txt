set(EXAMPLES
    bhatti_example_1dot4
)

foreach(var ${EXAMPLES})
    add_executable(${var} "${var}.cpp")
    target_compile_definitions(${var} PUBLIC USE_MKL)
    target_link_libraries(${var} PUBLIC MKL::MKL ${LACLIB_LIBS} fem2d)

    # this will run the example after build
    add_custom_target(run_${var} ALL ${var})
endforeach(var)
