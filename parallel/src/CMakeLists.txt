set(SOURCES main.c jacobi_parallel.c jacobi.h)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

add_executable(${CMAKE_PROJECT_NAME}_parallel ${SOURCES})
target_link_libraries(${CMAKE_PROJECT_NAME}_parallel m ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(${CMAKE_PROJECT_NAME}_parallel PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(${CMAKE_PROJECT_NAME}_parallel PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

