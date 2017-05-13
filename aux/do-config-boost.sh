#!/bin/bash -ex

# modify these paths for your system
BOOST_SRC_DIR=
BOOST_INSTALL_DIR=

cd ${BOOST_SRC_DIR}
./bootstrap.sh \
--with-libraries=signals,regex,filesystem,system,mpi,serialization,thread,program_options,exception \
--prefix=${BOOST_INSTALL_DIR} \
2>&1 | tee config_log
