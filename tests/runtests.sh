#!/bin/bash

set -e

echo --------------------------
echo units_test
./units_test

if test -f "./units20_test"; then
  echo --------------------------
  echo units20_test
  ./units20_test
fi

echo --------------------------
echo units_coords_test
./units_coords_test

if test -f "./eigen_units_test"; then
    echo --------------------------
    echo eigen_units_test
    ./eigen_units_test
fi
