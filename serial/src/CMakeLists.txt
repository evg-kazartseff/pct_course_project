set(SOURCES main.c jacobi.c jacobi.h)

add_executable(${CMAKE_PROJECT_NAME} ${SOURCES})
target_link_libraries(${CMAKE_PROJECT_NAME} m)