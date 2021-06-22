#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2021 Department of Computer Science, ETH Zurich
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



OBJS = opt_zones_closure.o opt_zones_incr_closure.o opt_zones_representation.o opt_mat.o opt_zones_predicate.o opt_zones_nary.o opt_zones_resize.o opt_zones_transfer.o

ifeq ($(IS_APRON),)
LIBS = -L../partitions_api -lpartitions -L../elina_auxiliary -lelinaux -L../elina_linearize -lelinalinearize $(MPFR_LIB_FLAG) -lmpfr $(GMP_LIB_FLAG) -lgmp -lm
INCLUDES = -I../ -I../elina_auxiliary -I../elina_linearize -I../partitions_api $(MPFR_INCLUDE_FLAG) $(GMP_INCLUDE_FLAG)
else
LIBS = -L../partitions_api -lpartitions -L$(APRON_PREFIX)/lib -lapron  -L../elina_linearize -lelinalinearize $(MPFR_LIB_FLAG) -lmpfr $(GMP_LIB_FLAG) -lgmp -lm
INCLUDES = -I$(APRON_PREFIX)/include -I../ -I../apron_interface -I../elina_linearize -I../partitions_api $(MPFR_INCLUDE_FLAG) $(GMP_INCLUDE_FLAG)
endif

SOINST = liboptzones.so

OPTZONESH = opt_zones.h opt_zones_internal.h opt_mat.h

all : liboptzones.so

opt_mat.o : opt_mat.h opt_mat.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_mat.o opt_mat.c $(LIBS)

opt_zones_predicate.o : opt_zones_predicate.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_zones_predicate.o opt_zones_predicate.c $(LIBS)

opt_zones_nary.o : opt_zones_nary.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_zones_nary.o opt_zones_nary.c $(LIBS)

opt_zones_resize.o : opt_zones_resize.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_zones_resize.o opt_zones_resize.c $(LIBS)

opt_zones_transfer.o : opt_zones_transfer.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_zones_transfer.o opt_zones_transfer.c $(LIBS)

opt_zones_closure.o : opt_zones_closure.h opt_zones_closure.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_zones_closure.o opt_zones_closure.c $(LIBS)

opt_zones_incr_closure.o : opt_zones_incr_closure.h opt_zones_incr_closure.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_zones_incr_closure.o opt_zones_incr_closure.c $(LIBS)

opt_zones_representation.o : opt_zones_representation.c opt_zones_internal.h
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_zones_representation.o opt_zones_representation.c $(LIBS)

liboptzones.so : $(OBJS) $(OPTZONESH)
	$(CC) -shared $(CC_ELINA_DYLIB) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o $(SOINST) $(OBJS) $(LIBS)

install:
	$(INSTALLd) $(LIBDIR); \
	for i in $(SOINST); do \
		$(INSTALL) $$i $(LIBDIR); \
	done; 

	$(INSTALLd) $(INCLDIR); \
	for i in $(OPTZONESH); do \
		$(INSTALL) $$i $(INCLDIR); \
	done; 

clean:
	-rm $(SOINST) 
	-rm *.o
