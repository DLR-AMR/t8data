FROM ubuntu:22.04
RUN apt-get update
# Install git, mpich, vim and zlib
RUN DEBIAN_FRONTEND=noninteractive apt-get install git mpich libz-dev vim -y
# Install t8code dependencies
RUN apt-get install libtool make automake -y
RUN mkdir source
RUN cd source && git clone https://holke@github.com/holke/t8code
RUN cd source/t8code && git checkout develop
RUN cd source/t8code && git submodule update --init && ./bootstrap
# Make folder sc, p4est, t8code_debug
RUN mkdir exec && mkdir exec/sc && mkdir exec/p4est && mkdir exec/t8code_debug 
# Make folder sc, p4est, t8code (release)
RUN mkdir exec/sc_release && mkdir exec/p4est_release && mkdir exec/t8code_release
# configure sc
RUN cd exec/sc && ../../source/t8code/sc/configure --enable-debug --enable-mpi --prefix=$PWD/install
# build sc
RUN cd exec/sc && make -j V=0 && make -j V=0 install
# configure sc release
RUN cd exec/sc_release && ../../source/t8code/sc/configure --enable-mpi --prefix=$PWD/install
# build sc release
RUN cd exec/sc_release && make -j V=0 && make -j V=0 install
# configure p4est
RUN cd exec/p4est && ../../source/t8code/p4est/configure --enable-debug --enable-mpi --prefix=$PWD/install --with-sc=$PWD/../sc/install
# build p4est
RUN cd exec/p4est && make -j V=0 && make -j V=0 install
# configure p4est release
RUN cd exec/p4est_release && ../../source/t8code/p4est/configure --enable-mpi --prefix=$PWD/install --with-sc=$PWD/../sc_release/install
# build p4est release
RUN cd exec/p4est_release && make -j V=0 && make -j V=0 install
# configure t8code
RUN cd exec/t8code_debug && ../../source/t8code/configure --enable-debug \
    --enable-mpi \
    CFLAGS="-Wall -pedantic -O3" \
    CXXFLAGS="-Wall -pedantic -O3" \
    --prefix=$PWD/install \
    --with-sc=$PWD/../sc/install \
    --with-p4est=$PWD/../p4est/install \
    --enable-static --disable-shared
# build t8code
RUN cd exec/t8code_debug && make clean
RUN cd exec/t8code_debug && make -j V=0 
RUN cd exec/t8code_debug make check -j V=0 && make install
# configure t8code release
RUN cd exec/t8code_release && ../../source/t8code/configure \
    --enable-mpi \
    CFLAGS="-Wall -pedantic -O3" \
    CXXFLAGS="-Wall -pedantic -O3" \
    --prefix=$PWD/install \
    --with-sc=$PWD/../sc_release/install \
    --with-p4est=$PWD/../p4est_release/install \
    --enable-static --disable-shared
# build t8code
RUN cd exec/t8code_release && make clean
RUN cd exec/t8code_release && make -j V=0
RUN cd exec/t8code_release make check -j V=0 && make install
