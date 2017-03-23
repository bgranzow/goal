### Introduction

To build Goal and its dependencies, you will need:
- [git][0]
- [CMake][1]
- A C++11 compiler with MPI wrappers (e.g. [MPICH3][2])

The Goal library has required dependencies on:
- [PUMI][3]
- [Trilinos][4]

### Building Trilinos

Troy makes use of the following Trilinos packages:
- [Teuchos][5] - for input file parsing and smart pointers
- [Pamgen][6] - for string parsing of math expressions
- [Sacado][7] - for automatic differentiation
- [Phalanx][8] - for local field evaluations
- [Belos][9] - for iterative linear algebra solvers
- [Ifpack2][10] - for ILU preconditioning
- [Zoltan][11] -for dynamic load balancing
- [MiniTensor][12] - for tensor algebra

To build these Trilinos packages, two external
dependencies must first be built:
- [boost][13]
- [parmetis][14]
- [yaml-cpp][15]

#### Building boost

```bash
wget https://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.tar.gz/download
tar -xzvf boost_1_61_0.tar.gz
cd boost_1_61_0
./bootstrap.sh --with-libraries=signals,regex,filesystem,system,mpi,serialization,thread,program_options,exception --prefix=<path-to-boost-install>
./b2 -j 4
./b2 install -j 4
```

Note: if your C++ mpi wrapper is not found as `mpicxx`, it may
be necessary to perform the step:

```bash
echo "using mpi : <path-to-mpicxx> ; " >> ~/user-config.jam
```

#### Building Parmetis

```bash
git clone https://github.com/ibaned/parmetis.git
cd parmetis
mkdir build
cd build
cmake .. \
-DCMAKE_INSTALL_PREFIX=<parmetis-install-dir> \
-DCMAKE_C_COMPILER=<mpi-c-wrapper> \
-DCMAKE_CXX_COMPILER=<mpi-c++-wrapper>
make -j 4
make install
```

#### Building yaml-cpp

```bash
git clone https://github.com/jbeder/yaml-cpp.git
cd yaml-cpp
mkdir build
cd build
cmake .. \
-DCMAKE_INSTALL_PREFIX=<yaml-cpp-install-dir>
make -j 4
make install
```

#### Building Trilinos

An example Trilinos [configuration script][16] is
provided with this source code. The first few lines
of this configuration script need to be modified
accordingly. To build Trilinos follow the steps:

```bash
git clone https://github.com/trilinos/Trilinos
cd Trilinos
mkdir build
cd build
* copy and modify do-config-trilinos.sh in this location *
./do-config-trilinos.sh
make -j 4
make install
```

### Building PUMI

An example PUMI [configuration script][17] is
provided with this source code. Additional information
about PUMI can be found [here][18]. to build PUMI
follow the steps:

```bash
git clone https://github.com/SCOREC/core.git
wget https://www.scorec.rpi.edu/pumi/pumi_test_meshes.tar.gz
tar xzf pumi_test_meshes.tar.gz
mkdir build
cd build
* copy and modify do-config-scorec.sh in this location *
./do-config-scorec.sh
make-j 4
make install
```

Note: the Zoltan dependency is provided by the
Trilinos installation, so the CMake Zoltan variables
in the PUMI configuration script should point to the
Trilinos install directory.

### Building Goal

Once the required dependencies have been installed,
Goal can be built by following the steps:

```bash
git clone https://bgranzow/goal
cd goal
mkdir build
cd build
cmake .. \
-DCMAKE_CXX_COMPILER=<path-to-mpicxx> \
-DCMAKE_INSTALL_PREFIX=<path-to-install>
make -j 4
make install
```

[0]:https://git-scm.com
[1]:https://cmake.org
[2]:https://www.mpich.org
[3]:https://github.com/scorec/core
[4]:https://github.com/trilinos/Trilinos
[5]:https://trilinos.org/packages/teuchos
[6]:https://trilinos.org/packages/pamgen
[7]:https://trilinos.org/packages/sacado
[8]:https://trilinos.org/packages/phalanx
[9]:https://trilinos.org/packages/belos
[10]:https://trilinos.org/packages/ifpack2
[11]:https://trilinos.org/packages/zoltan
[12]:https://github.com/trilinos/Trilinos/tree/master/packages/minitensor
[13]:http://www.boost.org
[14]:http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview
[15]:https://github.com/jbeder/yaml-cpp
[16]:https://github.com/bgranzow/goal/blob/master/aux/do-config-trilinos.sh
[17]:https://github.com/bgranzow/goal/blob/master/aux/do-config-scorec.sh
[18]:https://github.com/SCOREC/core/wiki/General-Build-instructions
