#!/bin/bash -ex

# Modify these paths for your system
SCOREC_DIR=/lore/granzb/core/install
TRILINOS_DIR=/lore/granzb/trilinos/install
GOAL_SRC_DIR=/lore/granzb/goal
GOAL_INSTALL_DIR=/lore/granzb/goal/install
GOAL_DATA_DIR=/lore/granzb/goal-data

cmake $GOAL_SRC_DIR \
-DCMAKE_CXX_COMPILER="mpicxx" \
-DCMAKE_INSTALL_PREFIX=$GOAL_INSTALL_DIR \
-DBUILD_TESTING=ON \
-DScorec_PREFIX=$SCOREC_DIR \
-DTrilinos_PREFIX=$TRILINOS_DIR \
-DGoal_DATA=$GOAL_DATA_DIR \
-DGoal_FAD_SIZE=40 \
-DGoal_ENABLE_ALL_APPS=ON \
-DGoal_EXTRA_CXX_FLAGS="-Wno-unused-parameter" \
2>&1 | tee config_log
