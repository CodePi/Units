default : units_test units_coords

units_test : units_test.cpp ../units.h
	$(CXX) -o $@ $< -std=c++11

units_coords_test : ../units_coords.h ../units.h
	$(CXX) -o $@ $< -std=c++11

# optional to build eigen_units_test
EIGEN_INC_DIR=/usr/include/eigen3
eigen_units_test : eigen_units_test.cpp ../units.h ../eigen_units.h
	$(CXX) -o $@ $< -std=c++11 -I$(EIGEN_INC_DIR)

all : default eigen_units_test

clean :
	rm -f units_test units_coords eigen_units_test
