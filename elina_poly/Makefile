#	
#	Copyright 2016 Software Reliability Lab, ETH Zurich	
#	Licensed under the Apache License, Version 2.0 (the "License"); 
#	you may not use this file except in compliance with the License.
#	You may obtain a copy of the License at
#
#		http://www.apache.org/licenses/LICENSE-2.0
#
#	Unless required by applicable law or agreed to in writing, software
#	distributed under the License is distributed on an "AS IS" BASIS,
#	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#	See the License for the specific language governing permissions and
#	limitations under the License.
#


include ../Makefile.config


AFLAGS := -D_GNU_SOURCE -pthread -fno-tree-vectorize -m64 -march=native

DFLAGS := -g  -DNUM_LONGLONGRAT -DTIMING -DHAS_APRON $(AFLAGS)

PREFIX = $(APRON_PREFIX)

LIBDIR = $(PREFIX)/lib

INCLDIR = $(PREFIX)/include 

INSTALL = install
INSTALLd = install -d

OBJS = opt_mf_qsort.o opt_pk_internal.o opt_pk_matrix.o opt_pk_test.o opt_pk_vector.o opt_pk_representation.o opt_pk_project.o opt_pk_constructor.o opt_pk_meetjoin.o opt_pk_assign.o opt_pk_widening.o opt_pk_resize.o opt_pk_expandfold.o opt_pk_extract.o opt_pk_user.o opt_pk_bit.o opt_pk_satmat.o opt_pk_cherni.o

INCLUDES = \
-I../newpolka \
-I../apron \
-I../mlapronidl \
-I../num \
-I../itv \
-I../octagons \
-I$(MPFR_PREFIX)/include -I$(GMP_PREFIX)/include \
-I../partitions_api

LIBS = -L../partitions_api -lpartitions -L../apron -L$(APRON_PREFIX) -lapron -L$(MPFR_PREFIX)/lib -lmpfr -L$(GMP_PREFIX)/lib -lgmp -lm

INSTALL = install
INSTALLd = install -d

SOINST = liboptpoly.so

OPTPOLYH = opt_mf_qsort.h opt_pk_matrix.h  opt_pk_internal.h opt_pk_test.h opt_pk_config.h opt_pk.h apron_wrapper.h

all : liboptpoly.so

opt_mf_qsort.o : opt_mf_qsort.h opt_mf_qsort.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_mf_qsort.o opt_mf_qsort.c $(LIBS)

opt_pk_internal.o : opt_pk_config.h opt_pk_internal.h opt_pk_internal.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_internal.o opt_pk_internal.c $(LIBS)

opt_pk_vector.o : opt_pk_vector.h opt_pk_vector.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_vector.o opt_pk_vector.c $(LIBS)	

opt_pk_matrix.o : opt_pk_matrix.h opt_pk_matrix.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_matrix.o opt_pk_matrix.c $(LIBS)

opt_pk_representation.o : opt_pk_representation.h opt_pk_representation.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_representation.o opt_pk_representation.c $(LIBS)

opt_pk_constructor.o : opt_pk_constructor.h opt_pk_constructor.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_constructor.o opt_pk_constructor.c $(LIBS)

opt_pk_project.o : opt_pk_project.c opt_pk_project.h
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_project.o opt_pk_project.c $(LIBS)

opt_pk_meetjoin.o : opt_pk_meetjoin.h opt_pk_representation.h opt_pk_meetjoin.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_meetjoin.o opt_pk_meetjoin.c $(LIBS)

opt_pk_assign.o : opt_pk_assign.h opt_pk_assign.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_assign.o opt_pk_assign.c $(LIBS)

opt_pk_widening.o : opt_pk_widening.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_widening.o opt_pk_widening.c $(LIBS)

opt_pk_expandfold.o : opt_pk_expandfold.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_expandfold.o opt_pk_expandfold.c $(LIBS)

opt_pk_extract.o : opt_pk_extract.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_extract.o opt_pk_extract.c $(LIBS)

opt_pk_resize.o : opt_pk_resize.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_resize.o opt_pk_resize.c $(LIBS)

opt_pk_test.o : opt_pk_test.h opt_pk_test.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_test.o opt_pk_test.c $(LIBS)

opt_pk_user.o : opt_pk_user.h opt_pk_user.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_user.o opt_pk_user.c $(LIBS)

opt_pk_bit.o : opt_pk_bit.h opt_pk_bit.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_bit.o opt_pk_bit.c $(LIBS)

opt_pk_satmat.o : opt_pk_satmat.h opt_pk_satmat.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_satmat.o opt_pk_satmat.c $(LIBS)

opt_pk_cherni.o : opt_pk_cherni.h opt_pk_cherni.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o opt_pk_cherni.o opt_pk_cherni.c $(LIBS)

liboptpoly.so : $(OBJS) $(OPTPOLYH)
	$(CC) -shared $(CFLAGS) $(DFLAGS) $(INCLUDES) -o $(SOINST) $(OBJS) $(LIBS)

install:
	$(INSTALLd) $(LIBDIR); \
	for i in $(SOINST); do \
		$(INSTALL) $$i $(LIBDIR); \
	done;

	$(INSTALLd) $(INCLDIR); \
	for i in $(OPTPOLYH); do \
		$(INSTALL) $$i $(INCLDIR); \
	done; 


clean:
	-rm *.o
	-rm *.so
