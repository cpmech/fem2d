set(EXAMPLES
    bhatti_five_member_truss_1dot4
    bhatti_plane_stress_1dot6
    felippa_three_member_truss
    smith_plane_strain_5dot2
)

foreach(var ${EXAMPLES})
    add_executable(${var} "${var}.cpp")
    target_compile_definitions(${var} PUBLIC USE_MKL)
    target_link_libraries(${var} PUBLIC MKL::MKL ${LACLIB_LIBS} fem2d)

    # this will run the example after build
    add_custom_target(run_${var} ALL ${var})
endforeach(var)
