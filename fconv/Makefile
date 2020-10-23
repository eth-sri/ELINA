#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
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

LIBDIR = $(ELINA_PREFIX)/lib
INCLDIR = $(ELINA_PREFIX)/include

INSTALL = install
INSTALLd = install -d

OBJS = fconv.o relaxation.o decomposition.o split_in_quadrants.o octahedron.o quadrants.o pdd.o sparse_cover.o mpq.o utils.o dynamic_bitset.o fp_mat.o curve_bounds.o

LIBS = $(GMP_LIB_FLAG) -lgmp -lcddgmp
INCLUDES = $(GMP_INCLUDE_FLAG)

ifneq ($(CDD_PREFIX),)
	INCLUDES += -I$(CDD_PREFIX)
endif

ifneq ($(BOOST_PREFIX),)
	INCLUDES += -I$(BOOST_PREFIX)
endif

SOINST = libfconv.so
FCONVH = fconv.h

CXXFLAGS := -std=c++14 -Wall -Wextra -Wpedantic -flto -mfma -mavx2 -DGMPRATIONAL -fPIC
ifeq ($(DEBUG),)
	CXXFLAGS := $(CXXFLAGS) -DNDEBUG -O3 -march=native
endif

all: libfconv.so tests

$(SOINST): $(OBJS)
	$(CXX) -shared $(CC_ELINA_DYLIB) $(CXXFLAGS) $(INCLUDES) -o $(SOINST) $(OBJS) $(LIBS)

tests: tests.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

curve_bounds.o : ../fppoly/curve_bounds.h ../fppoly/curve_bounds.c
	$(CXX) -c $(CXXFLAGS) -o curve_bounds.o ../fppoly/curve_bounds.c

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CDDFLAGS) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

install:
	$(INSTALLd) $(LIBDIR); \
	for i in $(SOINST); do \
		$(INSTALL) $$i $(LIBDIR); \
	done;

	$(INSTALLd) $(INCLDIR); \
	for i in $(FCONVH); do \
		$(INSTALL) $$i $(INCLDIR); \
	done;

clean:
	-rm *.o
	-rm *.so
	-rm tests

.PHONY: test
test: tests
	./tests
