#!/bin/bash -ex

# Modify these paths for your system
SCOREC_DIR=/home/bng/codes/core/install
TRILINOS_DIR=/home/bng/codes/trilinos/install
GOAL_SRC_DIR=/home/bng/codes/goal
GOAL_INSTALL_DIR=/home/bng/codes/goal/install
GOAL_DATA_DIR=/home/bng/codes/goal-data

cmake $GOAL_SRC_DIR \
-DCMAKE_CXX_COMPILER="mpicxx" \
-DCMAKE_INSTALL_PREFIX=$GOAL_INSTALL_DIR \
-DBUILD_TESTING=ON \
-DSCOREC_PREFIX=$SCOREC_DIR \
-DTrilinos_PREFIX=$TRILINOS_DIR \
-DGoal_DATA=$GOAL_DATA_DIR \
-DGoal_FAD_SIZE=40 \
2>&1 | tee config_log
