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

PREFIX = $(APRON_PREFIX)
LIBDIR = $(PREFIX)/lib

export CLASSPATH=.:/usr/local/apron.jar:/usr/local/gmp.jar:$CLASSPATH
############################

IFLAGS = \
$(GMP_INCLUDE_FLAG) \
$(MPFR_INCLUDE_FLAG) \
-I$(APRON_SRCROOT)/apron \
-I$(APRON_SRCROOT)/japron \
-I$(APRON_SRCROOT)/japron/apron \
-I$(APRON_SRCROOT)/japron/gmp \
-I../elina_oct \
-I../elina_poly \
-I../partitions_api \
-I$(JAVA_HOME)/include \
-I$(JAVA_HOME)/include/linux \
$(JNIINC) 

LFLAGS = \
$(GMP_LIB_FLAG) \
$(MPFR_LIB_FLAG) \
-L$(APRON_SRCROOT)/apron \
-L$(APRON_SRCROOT)/japron \
-L$(APRON_SRCROOT)/japron/apron \
-L$(APRON_SRCROOT)/japron/gmp \
-L. \
-L../elina_oct \
-L../elina_poly 


CFLAGSN = $(CFLAGS) -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -Wno-bad-function-cast

############################
#SOURCES
########

GMPLIBS   = -lmpfr -lgmp -lm


#ELINALIBS = -ljgmp -lboxD -loctD -lpolkaMPQ -lpolkaRll  -lapron -loptoct -loptpoly $(GMPLIBS)
ELINALIBS =  -ljapron -lapron -loptoct -loptpoly  $(GMPLIBS)

# .java / .c pairs
ELINAMODS =  OptPoly OptOctagon

ELINAJAVA = Test

############################
#OBJECTS
########

ELINAJ     = $(addprefix elina/,$(addsuffix .java,$(ELINAMODS) $(ELINAJAVA)))
ELINACLASS = $(subst .java,.class,$(ELINAJ))
ELINAC     = $(addprefix elina/elina_,$(addsuffix .c,$(ELINAMODS)))
ELINAO     = $(addprefix elina/elina_,$(addsuffix .o,$(ELINAMODS)))
ELINAH     = $(addprefix elina/elina_,$(addsuffix .h,$(ELINAMODS)))

ELINALL   = $(ELINACLASS) $(ELINAH) $(ELINAO) libjelina.so


# built and installed in JAVA_PREFIX
JAVAINST   = elina.jar 

APRONJAR   = apron.jar

# built and installed in LIBDIR
SOINST     = libjelina.so

############################
#TARGETS
########

all: $(ELINALL) $(JAVAINST) $(SOINST)


$(ELINACLASS): $(ELINAJ)
	$(JAVAC)  $+ 

elina/%.o:    elina/%.c $(ELINAH) 
	$(CC) $(CFLAGSN) -c $(IFLAGS) -Igmp $< -o $@

libjelina.so: $(ELINAO)
	$(CC) $(CFLAGSN) -shared $+ -o $@ $(LFLAGS) $(ELINALIBS)

elina/elina_%.h: elina/%.java
	$(JAVAH) -o $@ elina.$* 
	touch $@

elina.jar: $(ELINACLASS)
	$(JAR) cf $@ $+
	$(JAR) i $@


install:
	$(INSTALLd) $(LIBDIR); \
	for i in $(SOINST); do \
		$(INSTALL) $$i $(LIBDIR); \
	done; \
	$(INSTALLd) $(JAVA_PREFIX); \
	for i in $(JAVAINST); do \
		$(INSTALL) $$i $(JAVA_PREFIX); \
	done

uninstall:
	for i in $(SOINST); do \
		/bin/rm -f $(LIBDIR)/$$i; \
	done; \
	for i in $(JAVAINST); do \
		/bin/rm -f $(JAVA_PREFIX)/$$i; \
	done

clean:
	rm -f *~ \#* */*~ *.o *.class *.so *.jar */\#* */*.o */*.class */*.so  $(ELINAH)
	rm -rf doc


