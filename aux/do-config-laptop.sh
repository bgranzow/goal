# Modify these paths for your system
SCOREC_DIR=/home/bng/codes/core/install
TRILINOS_DIR=/home/bng/codes/trilinos/install
GOAL_INSTALL=/home/bng/codes/goal/install
GOAL_DATA_DIR=/home/bng/codes/goal-data

cmake \
 -D CMAKE_CXX_COMPILER="mpicxx" \
 -D CMAKE_INSTALL_PREFIX=$GOAL_INSTALL \
 -D BUILD_TESTING=ON \
 -D Scorec_PREFIX=$SCOREC_DIR \
 -D Trilinos_PREFIX=$TRILINOS_DIR \
 -D Goal_DATA=$GOAL_DATA_DIR \
 -D Goal_FAD_SIZE=40 \
 -D Goal_ENABLE_ALL_APPS=ON \
..
