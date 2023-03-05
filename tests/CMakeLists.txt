project(Units)

set(CMAKE_CXX_STANDARD 11)

# build units test
add_executable(units_test units_test.cpp)
add_test(units_test units_test)

# build units20 test
if(CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 11)
    add_executable(units20_test units20_test.cpp)
    set_property(TARGET units20_test PROPERTY CXX_STANDARD 20)
    add_test(units20_test units20_test)
endif()

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