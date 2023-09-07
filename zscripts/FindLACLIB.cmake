set(LACLIB_INCLUDE_SEARCH_PATH
    /usr/local/include/laclib
)

set(LACLIB_LIBRARY_SEARCH_PATH
    /usr/local/lib/laclib
)

find_path(LACLIB_INC laclib.h ${LACLIB_INCLUDE_SEARCH_PATH} NO_DEFAULT_PATH)
find_library(LACLIB_LIB NAMES laclib_intel PATHS ${LACLIB_LIBRARY_SEARCH_PATH} NO_DEFAULT_PATH)

set(LACLIB_FOUND 1)
foreach(var LACLIB_INC LACLIB_LIB)
    if(NOT ${var})
        set(LACLIB_FOUND 0)
    endif(NOT ${var})
endforeach(var)

if(LACLIB_FOUND)
    set(LACLIB_INCS ${LACLIB_INC})
    set(LACLIB_LIBS ${LACLIB_LIB})
else()
    message(FATAL_ERROR "cannot find laclib")
endif()
