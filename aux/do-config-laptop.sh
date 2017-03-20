# Modify these paths for your system
SCOREC_DIR=/Users/bng/codes/core/install-opt
TRILINOS_DIR=/Users/bng/codes/trilinos/install-opt
GOAL_INSTALL=/Users/bng/codes/goal/install-opt

cmake \
 -D CMAKE_CXX_COMPILER="mpicxx" \
 -D CMAKE_INSTALL_PREFIX=$GOAL_INSTALL \
 -D BUILD_TESTING=ON \
 -D Scorec_PREFIX=$SCOREC_DIR \
 -D Trilinos_PREFIX=$TRILINOS_DIR \
 -D GOAL_FAD_SIZE=40 \
 -D Goal_CXX_FLAGS="-std=c++11 -g -O2 -Werror -Wall -Wextra -ferror-limit=2" \
..
