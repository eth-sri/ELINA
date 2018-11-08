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

LIBDIR = $(ELINA_PREFIX)/lib

INCLDIR = $(ELINA_PREFIX)/include 

INSTALL = install
INSTALLd = install -d

OBJS = zonoml_internal.o zonoml_fun.o zonoml_reduced_product.o

ifeq ($(IS_APRON),)
LIBS = -L../partitions_api -lpartitions -L../elina_auxiliary -lelinaux -L../elina_linearize -lelinalinearize -L../elina_ml  -L../elina_zonotope -lzonotope $(MPFR_LIB_FLAG) -lmpfr $(GMP_LIB_FLAG) -lgmp -lm
INCLUDES = -I../ -I../elina_auxiliary -I../elina_linearize -I../elina_ml -I../elina_zonotope -I../partitions_api $(MPFR_INCLUDE_FLAG) $(GMP_INCLUDE_FLAG)
else
LIBS = -L../partitions_api -lpartitions -L$(APRON_PREFIX)/lib -lapron  -L../elina_linearize -lelinalinearize -L../elina_ml  -L../elina_zonotope -lzonotope $(MPFR_LIB_FLAG) -lmpfr $(GMP_LIB_FLAG) -lgmp -lm
INCLUDES = -I../apron_interface -I../ -I$(APRON_PREFIX)/include -I../elina_linearize -I../elina_ml -I../elina_zonotope -I../partitions_api $(MPFR_INCLUDE_FLAG) $(GMP_INCLUDE_FLAG)
endif

INSTALL = install
INSTALLd = install -d

SOINST = libzonoml.so

zonomlH = zonoml.h 

all : libzonoml.so elina_test_zonoml

libzonoml.so : $(OBJS) $(zonomlH)
	$(CC) -shared $(CC_ELINA_DYLIB) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o $(SOINST) $(OBJS) $(LIBS)


zonoml_internal.o : zonoml_internal.h zonoml_internal.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o zonoml_internal.o zonoml_internal.c $(LIBS)

zonoml_fun.o : zonoml_fun.h zonoml_fun.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o zonoml_fun.o zonoml_fun.c $(LIBS)

zonoml_reduced_product.o : zonoml_reduced_product.h zonoml_reduced_product.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o zonoml_reduced_product.o zonoml_reduced_product.c $(LIBS)

elina_test_zonoml : elina_test_zonoml.c libzonoml.so
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_test_zonoml elina_test_zonoml.c $(LIBS) -L. -lzonoml  -lzonotope -lelinaux -lelinalinearize -lpartitions -lmpfr -lgmp -lm

install:
	$(INSTALLd) $(LIBDIR); \
	for i in $(SOINST); do \
		$(INSTALL) $$i $(LIBDIR); \
	done;

	$(INSTALLd) $(INCLDIR); \
	for i in $(zonomlH); do \
		$(INSTALL) $$i $(INCLDIR); \
	done; 


clean:
	-rm *.o
	-rm *.so
	-rm elina_test_zonoml

