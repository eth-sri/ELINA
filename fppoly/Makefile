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

INSTALL = install
INSTALLd = install -d

OBJS = fppoly.o backsubstitute.o compute_bounds.o expr.o relu_approx.o leakyrelu_approx.o round_approx.o clip_approx.o batch_normalization.o sign_approx.o s_curve_approx.o parabola_approx.o log_approx.o pool_approx.o lstm_approx.o maxpool_convex_hull.o

ifeq ($(IS_APRON),)
LIBS = -L../partitions_api -lpartitions -L../elina_auxiliary -lelinaux -L../elina_linearize -lelinalinearize  -L../elina_zonotope -lzonotope $(MPFR_LIB_FLAG) -lmpfr $(GMP_LIB_FLAG) -lgmp -lm -lpthread -lcdd
INCLUDES = -I../elina_auxiliary -I../elina_linearize -I../elina_zonotope -I../partitions_api -I/home/gagandeepsi/gurobi903/linux64/include/  -I$(CDD_PREFIX) $(MPFR_INCLUDE_FLAG) $(GMP_INCLUDE_FLAG)
else
LIBS = -L../partitions_api -lpartitions -L$(APRON_PREFIX)/lib -lapron  -L../elina_linearize -lelinalinearize  -L../elina_zonotope -lzonotope $(MPFR_LIB_FLAG) -lmpfr $(GMP_LIB_FLAG) -lgmp -lm -lpthread -lcdd
INCLUDES = -I../apron_interface -I$(APRON_PREFIX)/include -I../elina_linearize -I../elina_zonotope -I../partitions_api  -I$(CDD_PREFIX) $(MPFR_INCLUDE_FLAG) $(GMP_INCLUDE_FLAG)
endif

ifneq ($(USE_GUROBI),)
DFLAGS += -DGUROBI
LIBS += -lgurobi90
endif

INSTALL = install
INSTALLd = install -d

SOINST = libfppoly.so

FPPOLYH = fppoly.h 

all : libfppoly.so 

libfppoly.so : $(OBJS) $(FPPOLYH)
	$(CC) -shared $(CC_ELINA_DYLIB) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o $(SOINST) $(OBJS) $(LIBS)

fppoly.o : fppoly.h fppoly.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o fppoly.o fppoly.c $(LIBS)

backsubstitute.o : backsubstitute.h backsubstitute.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o backsubstitute.o backsubstitute.c $(LIBS)

compute_bounds.o : compute_bounds.h compute_bounds.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o compute_bounds.o compute_bounds.c $(LIBS)

expr.o : expr.h expr.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o expr.o expr.c $(LIBS)

relu_approx.o : relu_approx.h relu_approx.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o relu_approx.o relu_approx.c $(LIBS)

leakyrelu_approx.o : leakyrelu_approx.h leakyrelu_approx.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o leakyrelu_approx.o leakyrelu_approx.c $(LIBS)

round_approx.o : round_approx.h round_approx.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o round_approx.o round_approx.c $(LIBS)
		
clip_approx.o : clip_approx.h clip_approx.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o clip_approx.o clip_approx.c $(LIBS)

batch_normalization.o : batch_normalization.h batch_normalization.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o batch_normalization.o batch_normalization.c $(LIBS)

sign_approx.o : sign_approx.h sign_approx.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o sign_approx.o sign_approx.c $(LIBS)

s_curve_approx.o : s_curve_approx.h s_curve_approx.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o s_curve_approx.o s_curve_approx.c $(LIBS)

parabola_approx.o : parabola_approx.h parabola_approx.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o parabola_approx.o parabola_approx.c $(LIBS)

log_approx.o : log_approx.h log_approx.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o log_approx.o log_approx.c $(LIBS)

pool_approx.o : pool_approx.h pool_approx.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o pool_approx.o pool_approx.c $(LIBS)

lstm_approx.o : lstm_approx.h lstm_approx.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o lstm_approx.o lstm_approx.c $(LIBS)


maxpool_convex_hull.o : maxpool_convex_hull.h maxpool_convex_hull.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) -o maxpool_convex_hull.o maxpool_convex_hull.c $(LIBS)


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

