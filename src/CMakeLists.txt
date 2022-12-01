add_subdirectory(common)
add_subdirectory(cpu)
add_subdirectory(omp_cpu)

set(EWALD_SRC
    main.cpp)

find_package(GSL)

if(GSL_FOUND)
    message("Using GSL")
    include_directories(${GSL_INCLUDE_DIRS})
endif()

if (EWALD_CPU)
    add_executable(ewald ${EWALD_SRC})
    target_link_libraries(ewald ewald_common ewald_cpu)
endif()
