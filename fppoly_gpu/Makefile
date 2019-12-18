#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2018 Department of Computer Science, ETH Zurich
#  This software is distributed under GNU Lesser General Public License Version 3.0.
#  For more information, see the ELINA project website at:
#  http://elina.ethz.ch
#
#  THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
#  EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
#  THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
#  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
#  TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY     
#  DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
#  SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
#  ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
#  CONTRACT, TORT OR OTHERWISE).
#
#



include ../Makefile.config
CXXFLAGS = -Wcast-qual -Wswitch -Wall -Wextra -Wundef -Wcast-align -fPIC -O3 -DNDEBUG -Werror-implicit-function-declaration -std=c++14 -DNUM_LONGLONGRAT -DNUM_DOUBLE -DTHRESHOLD=0.75 -DTIMING -fno-tree-vectorize -m64 -march=native

LIBDIR = $(ELINA_PREFIX)/lib

INCLDIR = $(ELINA_PREFIX)/include 

INSTALL = install
INSTALLd = install -d

OBJS = fppoly_gpu.o

ifeq ($(IS_APRON),)
LIBS = -L../elina_auxiliary -lelinaux -L../elina_linearize -lelinalinearize
INCLUDES = -I../elina_auxiliary -I../elina_linearize
else
LIBS = -L$(APRON_PREFIX)/lib -lapron  -L../elina_linearize -lelinalinearize
INCLUDES = -I../apron_interface -I$(APRON_PREFIX)/include -I../elina_linearize
endif

INSTALL = install
INSTALLd = install -d

SOINST = libfppoly_gpu.so

FPPOLYH = fppoly_gpu.h

all : libfppoly_gpu.so

libfppoly_gpu.so : $(OBJS) $(FPPOLYH)
	nvcc -ccbin g++ -shared $(CC_ELINA_DYLIB) --compiler-options "$(CXXFLAGS)" -arch=sm_61 $(INCLUDES) -o $(SOINST) $(OBJS) $(LIBS)


fppoly_gpu.o : cuda_util.h rounding.h interval.h linear.h relu_approx.h fppoly_gpu.h fppoly_gpu.cu
	nvcc -ccbin g++ -c --compiler-options "$(CXXFLAGS)" -arch=sm_61 $(INCLUDES) -o fppoly_gpu.o fppoly_gpu.cu $(LIBS)



install:
	$(INSTALLd) $(LIBDIR); \
	for i in $(SOINST); do \
		$(INSTALL) $$i $(LIBDIR); \
	done;

	$(INSTALLd) $(INCLDIR); \
	for i in $(FPPOLYH); do \
		$(INSTALL) $$i $(INCLDIR); \
	done; 


clean:
	-rm *.o
	-rm *.so

