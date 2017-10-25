#!/bin/bash -ex

# Modify these paths for your system
SCOREC_DIR=/Users/bng/codes/core/install-opt
TRILINOS_DIR=/Users/bng/codes/trilinos/install-opt
GOAL_SRC_DIR=/Users/bng/codes/goal
GOAL_INSTALL_DIR=/Users/bng/goal/install-opt

cmake $GOAL_SRC_DIR \
-DCMAKE_CXX_COMPILER="mpicxx" \
-DCMAKE_INSTALL_PREFIX=$GOAL_INSTALL_DIR \
-DBUILD_TESTING=ON \
-DSCOREC_PREFIX=$SCOREC_DIR \
-DTrilinos_PREFIX=$TRILINOS_DIR \
-DGOAL_FAD_SIZE=16 \
-DGOAL_ENABLE_SNAPPING=OFF \
-DGOAL_EXTRA_CXX_FLAGS="-fno-omit-frame-pointer" \
2>&1 | tee config_log
