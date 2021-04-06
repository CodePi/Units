project(Units)

set(CMAKE_CXX_STANDARD 11)

# build units test
add_executable(units_test units_test.cpp)

# build eigen_units_test, if eigen is available
find_package (Eigen3)
string(COMPARE NOTEQUAL "${EIGEN3_INCLUDE_DIRS}" "" eigen_include_set)
if(${eigen_include_set})
    add_executable(eigen_units_test eigen_units_test.cpp)
    target_include_directories(eigen_units_test PUBLIC ${EIGEN3_INCLUDE_DIRS})
    add_test(eigen_units_test eigen_units_test)
endif()

add_executable(units_coords_test units_coords_test.cpp)
add_test(units_coords_test units_coords_test)