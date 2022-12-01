set(EWALD_CPU_SRC
    functions.cpp)

add_library(ewald_cpu ${EWALD_CPU_SRC})
target_link_libraries(ewald_cpu ewald_common)
