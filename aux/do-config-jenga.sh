#!/bin/bash -ex

# Modify these paths for your system
TRILINOS_DIR=/usr/local/trilinos/latest
GOAL_SRC_DIR=/lore/granzb/goal
GOAL_INSTALL_DIR=/lore/granzb/goal/install

cmake $GOAL_SRC_DIR \
-DCMAKE_CXX_COMPILER="mpicxx" \
-DCMAKE_INSTALL_PREFIX=$GOAL_INSTALL_DIR \
-DBUILD_TESTING=ON \
-DSCOREC_PREFIX=$TRILINOS_DIR \
-DTrilinos_PREFIX=$TRILINOS_DIR \
-DGoal_FAD_SIZE=70 \
2>&1 | tee config_log
