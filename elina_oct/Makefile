#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2017 Department of Computer Science, ETH Zurich
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

#PREFIX = /usr/local

LIBDIR = $(ELINA_PREFIX)/lib

INCLDIR = $(ELINA_PREFIX)/include
#SOBJS = $(LIBDIR)/liboptoct.so 

ifeq ($(IS_VECTOR),)
CLOSURE_OBJS = opt_oct_closure_comp_sparse.o  opt_oct_incr_closure_comp_sparse.o opt_oct_closure_dense_scalar.o opt_oct_incr_closure_dense_scalar.o
CLOSURE_C = opt_oct_closure_comp_sparse.c  opt_oct_incr_closure_comp_sparse.c opt_oct_closure_dense_scalar.c opt_oct_incr_closure_dense_scalar.c
CLOSUREH = opt_oct_closure_comp_sparse.h  opt_oct_incr_closure_comp_sparse.h opt_oct_closure_dense_scalar.h opt_oct_incr_closure_dense_scalar.h 
else
CLOSURE_OBJS = opt_oct_closure_comp_sparse.o  opt_oct_incr_closure_comp_sparse.o opt_oct_closure_dense.o opt_oct_incr_closure_dense.o
CLOSURE_C = opt_oct_closure_comp_sparse.c  opt_oct_incr_closure_comp_sparse.c opt_oct_closure_dense.c opt_oct_incr_closure_dense.c
CLOSUREH =  opt_oct_closure_comp_sparse.h  opt_oct_incr_closure_comp_sparse.h opt_oct_closure_dense.h opt_oct_incr_closure_dense.h vector_intrin.h
endif


OBJS = $(CLOSURE_OBJS) opt_oct_nary.o opt_oct_resize.o opt_oct_predicate.o opt_oct_representation.o opt_oct_transfer.o opt_oct_hmat.o

ifeq ($(IS_APRON),)
LIBS = -L../partitions_api -lpartitions -L../elina_auxiliary -lelinaux -L../elina_linearize -lelinalinearize $(MPFR_LIB_FLAG) -lmpfr $(GMP_LIB_FLAG) -lgmp -lm
INCLUDES = -I../elina_auxiliary -I../elina_linearize -I../partitions_api $(MPFR_INCLUDE_FLAG) $(GMP_INCLUDE_FLAG)
else
LIBS = -L../partitions_api -lpartitions -L$(APRON_PREFIX)/lib -lapron  -L../elina_linearize -lelinalinearize $(MPFR_LIB_FLAG) -lmpfr $(GMP_LIB_FLAG) -lgmp -lm
INCLUDES = -I$(APRON_PREFIX)/include -I../apron_interface -I../elina_linearize -I../partitions_api $(MPFR_INCLUDE_FLAG) $(GMP_INCLUDE_FLAG)
endif


INSTALL = install
INSTALLd = install -d

SOINST = liboptoct.so

OPTOCTH = opt_oct.h opt_oct_internal.h opt_oct_hmat.h $(CLOSUREH)

all : liboptoct.so elina_test_oct


opt_oct_closure_comp_sparse.o : opt_oct_closure_comp_sparse.h opt_oct_closure_comp_sparse.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_oct_closure_comp_sparse.o opt_oct_closure_comp_sparse.c $(LIBS)

opt_oct_incr_closure_comp_sparse.o : opt_oct_incr_closure_comp_sparse.h opt_oct_incr_closure_comp_sparse.c opt_oct_closure_comp_sparse.o
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_oct_incr_closure_comp_sparse.o opt_oct_incr_closure_comp_sparse.c $(LIBS) opt_oct_closure_comp_sparse.o

#ifneq ($(IS_VECTOR),)
opt_oct_closure_dense_scalar.o : opt_oct_closure_dense_scalar.h opt_oct_closure_dense_scalar.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_oct_closure_dense_scalar.o opt_oct_closure_dense_scalar.c $(LIBS)

opt_oct_incr_closure_dense_scalar.o : opt_oct_incr_closure_dense_scalar.h opt_oct_incr_closure_dense_scalar.c opt_oct_closure_dense_scalar.o
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_oct_incr_closure_dense_scalar.o opt_oct_incr_closure_dense_scalar.c $(LIBS) opt_oct_closure_dense_scalar.o


#else
opt_oct_closure_dense.o : opt_oct_closure_dense.h opt_oct_closure_dense.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_oct_closure_dense.o opt_oct_closure_dense.c $(LIBS)

opt_oct_incr_closure_dense.o : opt_oct_incr_closure_dense.h opt_oct_incr_closure_dense.c opt_oct_closure_dense.o
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_oct_incr_closure_dense.o opt_oct_incr_closure_dense.c $(LIBS) opt_oct_closure_dense.o

#endif

opt_oct_hmat.o : opt_oct_hmat.h opt_oct_hmat.c 
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_oct_hmat.o opt_oct_hmat.c $(LIBS) 

opt_oct_nary.o : opt_oct_nary.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_oct_nary.o opt_oct_nary.c $(LIBS) 

opt_oct_resize.o : opt_oct_resize.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_oct_resize.o opt_oct_resize.c $(LIBS) 

opt_oct_predicate.o : opt_oct_predicate.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_oct_predicate.o opt_oct_predicate.c $(LIBS) 

opt_oct_representation.o : opt_oct_representation.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_oct_representation.o opt_oct_representation.c $(LIBS) 

opt_oct_transfer.o : opt_oct_transfer.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_oct_transfer.o opt_oct_transfer.c $(LIBS) 

liboptoct.so : $(OBJS) $(OPTOCTH)
	$(CC) -shared -install_name $(LIBDIR)/$@ $(CFLAGS) $(DFLAGS) $(INCLUDES) -o $(SOINST) $(OBJS) $(LIBS)

elina_test_oct : elina_test_oct.c liboptoct.so
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_test_oct elina_test_oct.c $(LIBS) -L. -loptoct

install:
	$(INSTALLd) $(LIBDIR); \
	for i in $(SOINST); do \
		$(INSTALL) $$i $(LIBDIR); \
	done; 

	$(INSTALLd) $(INCLDIR); \
	for i in $(OPTOCTH); do \
		$(INSTALL) $$i $(INCLDIR); \
	done; 

clean:
	-rm $(SOINST) 
	-rm *.o
	-rm elina_test_oct
