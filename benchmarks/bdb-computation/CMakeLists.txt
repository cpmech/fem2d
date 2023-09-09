add_executable(bmark_bdb "main.cpp")
target_compile_definitions(bmark_bdb PUBLIC USE_MKL)
target_link_libraries(bmark_bdb PUBLIC MKL::MKL ${LACLIB_LIBS} fem2d)
