
srcs := voxelImageProcess.cpp   voxelToFoamPar.cpp  voxelToFoam.cpp
tsts := test.py
all: $(srcs)
$(info srcs: $(srcs),   tsts: $(tsts) )


msSrc ?= $(abspath ..)
USE_ZLIB=1
USE_TIFF=1
USE_OMP=1
USE_CPP17=1
USE_SINGLECPP=1
USE_msTEST=1
include  ../Makefile.in

