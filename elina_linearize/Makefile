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

AFLAGS := -D_GNU_SOURCE -pthread -fno-tree-vectorize -m64 -march=native

DFLAGS := -g $(AFLAGS)

LIBDIR = $(ELINA_PREFIX)/lib

INCLDIR = $(ELINA_PREFIX)/include 

INSTALL = install
INSTALLd = install -d

OBJS = elina_scalar_arith.o elina_interval_arith.o  elina_coeff_arith.o elina_linexpr0_arith.o elina_linearize.o elina_linearize_texpr.o

INCLUDES = -I$(MPFR_PREFIX)/include -I$(GMP_PREFIX)/include

LIBS = -L$(MPFR_PREFIX)/lib -lmpfr -L$(GMP_PREFIX)/lib -lgmp -lm  -lelinaux

INSTALL = install
INSTALLd = install -d

SOINST = libelinalinearize.so

ELINALINEARIZEH = elina_scalar_arith.h elina_interval_arith.h elina_coeff_arith.h elina_linexpr0_arith.h elina_linearize.h elina_linearize_texpr.h 

all :	libelinalinearize.so 

elina_scalar_arith.o : elina_scalar_arith.h elina_scalar_arith.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_scalar_arith.o elina_scalar_arith.c $(LIBS)

elina_interval_arith.o : elina_interval_arith.h elina_interval_arith.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_interval_arith.o elina_interval_arith.c $(LIBS)

elina_coeff_arith.o : elina_coeff_arith.h elina_coeff_arith.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_coeff_arith.o elina_coeff_arith.c $(LIBS)

elina_linexpr0_arith.o : elina_linexpr0_arith.h elina_linexpr0_arith.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_linexpr0_arith.o elina_linexpr0_arith.c $(LIBS)

elina_linearize.o : elina_linearize.h elina_linearize.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_linearize.o elina_linearize.c $(LIBS)

elina_linearize_texpr.o : elina_linearize_texpr.h elina_linearize_texpr.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_linearize_texpr.o elina_linearize_texpr.c $(LIBS)

libelinalinearize.so : $(OBJS) $(ELINALINEARIZEH)
	$(CC) -shared $(CFLAGS) $(DFLAGS) $(INCLUDES) -o $(SOINST) $(OBJS) $(LIBS)

install:
	$(INSTALLd) $(LIBDIR); \
	for i in $(SOINST); do \
		$(INSTALL) $$i $(LIBDIR); \
	done;

	$(INSTALLd) $(INCLDIR); \
	for i in $(ELINALINEARIZEH); do \
		$(INSTALL) $$i $(INCLDIR); \
	done; 


clean:
	-rm *.o
	-rm *.so
