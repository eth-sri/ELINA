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
include vars.mk

#---------------------------------------
# Programs
#---------------------------------------

#---------------------------------------
# Flags
#---------------------------------------

ICFLAGS += $(BASE_ICFLAGS) $(ML_ICFLAGS)
LDFLAGS += $(BASE_LIFLAGS)
CMXSINC = $(APRON_CMXSINC) -I .

APRON_CAML_INCLUDE = $(shell $(OCAMLFIND) query apron)


#---------------------------------------
# Files
#---------------------------------------

CAML_TO_INSTALL = \
elina_poly.idl elina_poly.mli elina_poly.ml elina_poly.cmi \
elina_poly.cma \
libelina_poly_caml.a \

ifneq ($(HAS_OCAMLOPT),)
CAML_TO_INSTALL += 							\
	$(call OCAMLOPT_TARGETS, $(addprefix elina_poly,bc))
endif
ifneq ($(HAS_SHARED),)
CAML_TO_INSTALL += dllelina_poly_caml.so
endif

LIBS = -lapron -lmpfr -lgmp -lm

#---------------------------------------
# Rules
#---------------------------------------

all: ml mlexample.byte

ml: elina_poly.mli elina_poly.ml elina_poly.cmi mllib 

mllib: elina_poly.cma libelina_poly_caml.a mlextra

ifneq ($(HAS_OCAMLOPT),)
mlextra: $(call OCAMLOPT_TARGETS, elina_polybc)
endif
ifneq ($(HAS_SHARED),)
mlextra: dllelina_poly_caml.so
endif


mlexample.byte: mlexample.ml elina_poly.cma
	$(OCAMLC) $(OCAMLFLAGS) -I $(MLGMPIDL_LIB) -I $(APRON_CAML_INCLUDE) -o $@ bigarray.cma gmp.cma apron.cma elina_poly.cma $<

clean:
	/bin/rm -f *.[ao] *.so
	/bin/rm -fr *.annot *.cm[ioax] *.cmx[as]
	/bin/rm -fr tmp
	/bin/rm -fr elina_poly_caml.* elina_poly.ml elina_poly.mli
	/bin/rm  mlexample.byte

distclean: clean

install:
	$(INSTALLd) $(APRON_INCLUDE) $(APRON_LIB)
	for i in $(CAML_TO_INSTALL); do \
		if test -f $$i; then $(INSTALL) $$i $(APRON_LIB); fi; \
	done

uninstall:
	/bin/rm -f $(APRON_LIB)/libelina_poly.a 
	/bin/rm -f $(APRON_LIB)/dllelina_poly.so
ifeq ($(OCAMLFIND),)
	/bin/rm -f $(APRON_LIB)/elina_poly.mli $(APRON_LIB)/elina_poly.ml $(APRON_LIB)/elina_poly.cm[ix] $(APRON_LIB)/elina_poly.idl $(APRON_LIB)/elina_poly*.cma $(APRON_LIB)/elina_poly*.cmx[as] $(APRON_LIB)/elina_poly*.a
endif


#---------------------------------------
# IMPLICIT RULES AND DEPENDENCIES
#---------------------------------------

.SUFFIXES: .tex .c .h .a .o

.PRECIOUS: libelina_poly_caml.a dllelina_poly_caml.so
.PRECIOUS: %.cmo %.cmx


#---------------------------------------
# OCaml binding
#---------------------------------------

dllelina_poly_caml.so libelina_poly_caml.a: elina_poly_caml.o
	$(OCAMLMKLIB) -o elina_poly_caml $< -L. -loptpoly $(LDFLAGS) $(LIBS)

#---------------------------------------
# OCaml rules
#---------------------------------------

elina_poly.cma: elina_poly.cmo libelina_poly_caml.a 
	$(OCAMLMKLIB) -o elina_poly$* -oc elina_poly$*_caml elina_poly.cmo -loptpoly $(LIBS)

elina_poly.cmxa elina_poly.a: elina_poly.cmx libelina_poly_caml.a 
	$(OCAMLMKLIB) -o elina_poly$* -oc elina_poly$*_caml elina_poly.cmx -loptpoly $(LIBS)


#---------------------------------------
# IDL rules
#---------------------------------------

%.ml %.mli %_caml.c: %.idl perlscript_caml.pl $(MLAPRONDIR)/manager.idl
	mkdir -p tmp
	cp $(MLAPRONDIR)/manager.idl $*.idl tmp/
	cd tmp && $(CAMLIDL) -no-include -nocpp $*.idl
	cp tmp/$*_stubs.c $*_caml.c
	$(PERL) perlscript_caml.pl < tmp/$*.ml >$*.ml
	$(PERL) perlscript_caml.pl < tmp/$*.mli >$*.mli

rebuild:
	@echo "make rebuild is no longer necessary"

#---------------------------------------
# C generic rules
#---------------------------------------

%.o: %.c
	$(CC) $(CFLAGS) $(ICFLAGS) -DNUMINT_MPQ -c -o $@ $<

#---------------------------------------
# ML generic rules
#---------------------------------------

%.cmi: %.mli
	$(OCAMLC) $(OCAMLFLAGS) $(OCAMLINC) -c $<

%.cmo: %.ml %.cmi
	$(OCAMLC) $(OCAMLFLAGS) $(OCAMLINC) -c $<

%.cmx: %.ml %.cmi
	$(OCAMLOPT) $(OCAMLOPTFLAGS) $(OCAMLINC) -c $<

%.cmxs: %.cmxa
	$(OCAMLOPT_CMXS) $(CMXSINC) -o $@ $<

