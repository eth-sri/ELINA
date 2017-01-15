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

LIBDIR = $(ELINA_PREFIX)/lib

INCLDIR = $(ELINA_PREFIX)/include 

INSTALL = install
INSTALLd = install -d

OBJS = elina_scalar.o elina_interval.o elina_coeff.o elina_dimension.o elina_linexpr0.o elina_lincons0.o elina_manager.o elina_abstract0.o elina_texpr0.o elina_tcons0.o

INCLUDES = -I$(MPFR_PREFIX)/include -I$(GMP_PREFIX)/include -I../elina_linearize

LIBS = -lmpfr -lgmp -lm -L../elina_linearize

INSTALL = install
INSTALLd = install -d

SOINST = libelinaux.so

ELINAAUXH = elina_config.h elina_scalar.h elina_interval.h elina_coeff.h elina_dimension.h elina_linexpr0.h elina_lincons0.h elina_manager.h elina_abstract0.h elina_texpr0.h elina_tcons0.h

all :	libelinaux.so

elina_scalar.o : elina_scalar.h elina_scalar.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_scalar.o elina_scalar.c $(LIBS)

elina_interval.o : elina_interval.h elina_interval.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_interval.o elina_interval.c $(LIBS)

elina_coeff.o : elina_coeff.h elina_coeff.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_coeff.o elina_coeff.c $(LIBS)	

elina_dimension.o : elina_dimension.h elina_dimension.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_dimension.o elina_dimension.c $(LIBS)

elina_linexpr0.o : elina_linexpr0.h elina_linexpr0.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_linexpr0.o elina_linexpr0.c $(LIBS)	

elina_lincons0.o : elina_lincons0.h elina_lincons0.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_lincons0.o elina_lincons0.c $(LIBS)

elina_manager.o : elina_manager.h elina_manager.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_manager.o elina_manager.c $(LIBS)	

elina_abstract0.o : elina_abstract0.h elina_abstract0.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_abstract0.o elina_abstract0.c $(LIBS)	

elina_texpr0.o : elina_texpr0.h elina_texpr0.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_texpr0.o elina_texpr0.c $(LIBS)	

elina_tcons0.o : elina_tcons0.h elina_tcons0.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o elina_tcons0.o elina_tcons0.c $(LIBS)


libelinaux.so : $(OBJS) $(ELINAAUXH)
	$(CC) -shared $(CFLAGS) $(DFLAGS) $(INCLUDES) -o $(SOINST) $(OBJS) $(LIBS)

install:
	$(INSTALLd) $(LIBDIR); \
	for i in $(SOINST); do \
		$(INSTALL) $$i $(LIBDIR); \
	done;

	$(INSTALLd) $(INCLDIR); \
	for i in $(ELINAAUXH); do \
		$(INSTALL) $$i $(INCLDIR); \
	done; 


clean:
	-rm *.o
	-rm *.so
