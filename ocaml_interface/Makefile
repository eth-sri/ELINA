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
include vars.mk

#---------------------------------------
# Programs
#---------------------------------------

#---------------------------------------
# Flags
#---------------------------------------

# add the paths of dynamically loaded libraries
export LD_LIBRARY_PATH=.:../elina_poly:../elina_oct:

ICFLAGS += $(BASE_ICFLAGS) $(ML_ICFLAGS)
LDFLAGS += $(BASE_LIFLAGS)
CMXSINC = $(APRON_CMXSINC) -I .

ICFLAGS += -I ../partitions_api -I ../elina_linearize -I ../elina_auxiliary -I ../elina_poly -I ../elina_oct

LDFLAGS += -L../partitions_api -L../elina_linearize  -L../elina_auxiliary -L../elina_poly -L../elina_oct -L.

ifeq ($(USE_OPAM),)
	LDFLAGS +=  -R../partitions_api -R../elina_linearize  -R../elina_auxiliary -R../elina_poly  -cclib -R../elina_oct -cclib -R.
endif


APRON_CAML_INCLUDE = $(shell $(OCAMLFIND) query apron)


#---------------------------------------
# Files
#---------------------------------------

CAML_TO_INSTALL = \
elina_poly.idl elina_poly.mli elina_poly.ml elina_poly.cmi \
elina_poly.a elina_poly.cma elina_poly.cmx elina_poly.cmxa \
libelina_poly_caml.a dllelina_poly.so dllelina_poly_caml.so \
elina_oct.idl elina_oct.mli elina_oct.ml elina_oct.cmi \
elina_oct.a elina_oct.cma elina_oct.cmx elina_oct.cmxa \
libelina_oct_caml.a dllelina_oct.so dllelina_oct_caml.so \


ifneq ($(HAS_OCAMLOPT),)
CAML_TO_INSTALL += 							\
	$(call OCAMLOPT_TARGETS, $(addprefix elina_poly,bc))
CAML_TO_INSTALL += 							\
	$(call OCAMLOPT_TARGETS, $(addprefix elina_oct,bc))
endif

LIBS = -lapron -lmpfr -lgmp -lm

#---------------------------------------
# Rules
#---------------------------------------

all: ml 

example: mlexample_poly.byte mlexample_poly.opt mlexample_oct.byte mlexample_oct.opt

ml: elina_poly.mli elina_poly.ml elina_poly.cmi elina_oct.mli elina_oct.ml elina_oct.cmi mllib 

mllib: elina_poly.cmxa elina_poly.cma elina_poly.a libelina_poly_caml.a elina_oct.cmxa elina_oct.cma elina_oct.a libelina_oct_caml.a mlextra

ifneq ($(HAS_OCAMLOPT),)
mlextra: $(call OCAMLOPT_TARGETS, elina_polybc) $(call OCAMLOPT_TARGETS, elina_octbc) 
endif



mlexample_poly.byte: mlexample_poly.ml elina_poly.cma
	$(OCAMLC) $(OCAMLFLAGS) -I $(MLGMPIDL_LIB) -I $(APRON_CAML_INCLUDE) -o $@ bigarray.cma gmp.cma apron.cma elina_poly.cma $<

mlexample_poly.opt: mlexample_poly.ml elina_poly.cmxa
	$(OCAMLOPT) $(OCAMLFLAGS) -I $(MLGMPIDL_LIB) -I $(APRON_CAML_INCLUDE) -o $@ bigarray.cmxa gmp.cmxa apron.cmxa elina_poly.cmxa $<

mlexample_oct.byte: mlexample_oct.ml elina_oct.cma
	$(OCAMLC) $(OCAMLFLAGS) -I $(MLGMPIDL_LIB) -I $(APRON_CAML_INCLUDE) -o $@ bigarray.cma gmp.cma apron.cma elina_oct.cma $<

mlexample_oct.opt: mlexample_oct.ml elina_oct.cmxa
	$(OCAMLOPT) $(OCAMLFLAGS) -I $(MLGMPIDL_LIB) -I $(APRON_CAML_INCLUDE) -o $@ bigarray.cmxa gmp.cmxa apron.cmxa elina_oct.cmxa $<

clean:
	/bin/rm -f *.[ao] *.so
	/bin/rm -fr *.annot *.cm[ioax] *.cmx[as]
	/bin/rm -fr tmp
	/bin/rm -fr elina_poly_caml.* elina_poly.ml elina_poly.mli
	/bin/rm -fr elina_oct_caml.* elina_oct.ml elina_oct.mli
	/bin/rm -f mlexample_poly.byte
	/bin/rm -f mlexample_poly.opt
	/bin/rm -f mlexample_oct.byte
	/bin/rm -f mlexample_oct.opt

distclean: clean

install:
	$(INSTALLd) $(ELINA_INCLUDE) $(ELINA_LIB)
	for i in $(CAML_TO_INSTALL); do \
		if test -f $$i; then $(INSTALL) $$i $(ELINA_LIB); fi; \
	done

uninstall:
	/bin/rm -f $(ELINA_LIB)/libelina_poly.a
	/bin/rm -f $(ELINA_LIB)/dllelina_poly.so
	/bin/rm -f $(ELINA_LIB)/libelina_oct.a
	/bin/rm -f $(ELINA_LIB)/dllelina_oct.so
ifeq ($(OCAMLFIND),)
	/bin/rm -f $(ELINA_LIB)/elina_poly.mli $(ELINA_LIB)/elina_poly.ml $(ELINA_LIB)/elina_poly.cm[ix] $(ELINA_LIB)/elina_poly.idl $(ELINA_LIB)/elina_poly*.cma $(ELINA_LIB)/elina_poly*.cmx[as] $(ELINA_LIB)/elina_poly*.a
	/bin/rm -f $(ELINA_LIB)/elina_oct.mli $(ELINA_LIB)/elina_oct.ml $(ELINA_LIB)/elina_oct.cm[ix] $(ELINA_LIB)/elina_oct.idl $(ELINA_LIB)/elina_oct*.cma $(ELINA_LIB)/elina_oct*.cmx[as] $(ELINA_LIB)/elina_oct*.a
endif


#---------------------------------------
# IMPLICIT RULES AND DEPENDENCIES
#---------------------------------------

.SUFFIXES: .tex .c .h .a .o

.PRECIOUS: libelina_poly_caml.a dllelina_poly_caml.so libelina_oct_caml.a dllelina_oct_caml.so
.PRECIOUS: %.cmo %.cmx


#---------------------------------------
# OCaml binding
#---------------------------------------

dllelina_poly_caml.so libelina_poly_caml.a: elina_poly_caml.o
	$(OCAMLMKLIB) -dllpath $(ELINA_PREFIX)/lib -L$(ELINA_PREFIX)/lib -o elina_poly_caml $< -L. -loptpoly $(LDFLAGS) $(LIBS)
	
dllelina_oct_caml.so libelina_oct_caml.a: elina_oct_caml.o
	$(OCAMLMKLIB) -dllpath $(ELINA_PREFIX)/lib -L$(ELINA_PREFIX)/lib -o elina_oct_caml $< -L. -loptoct $(LDFLAGS) $(LIBS)

#---------------------------------------
# OCaml rules
#---------------------------------------

elina_poly.cma: elina_poly.cmo libelina_poly_caml.a 
	$(OCAMLMKLIB) -dllpath $(ELINA_PREFIX)/lib -L$(ELINA_PREFIX)/lib -o elina_poly$* -oc elina_poly$*_caml elina_poly.cmo -loptpoly $(LIBS)

elina_poly.cmxa elina_poly.a: elina_poly.cmx libelina_poly_caml.a 
	$(OCAMLMKLIB) -dllpath $(ELINA_PREFIX)/lib -L$(ELINA_PREFIX)/lib -o elina_poly$* -oc elina_poly$*_caml elina_poly.cmx -loptpoly $(LIBS)


elina_oct.cma: elina_oct.cmo libelina_oct_caml.a 
	$(OCAMLMKLIB) -dllpath $(ELINA_PREFIX)/lib -L$(ELINA_PREFIX)/lib -o elina_oct$* -oc elina_oct$*_caml elina_oct.cmo -loptoct $(LIBS)

elina_oct.cmxa elina_oct.a: elina_oct.cmx libelina_oct_caml.a 
	$(OCAMLMKLIB) -dllpath $(ELINA_PREFIX)/lib -L$(ELINA_PREFIX)/lib -o elina_oct$* -oc elina_oct$*_caml elina_oct.cmx -loptoct $(LIBS)

#---------------------------------------
# IDL rules
#---------------------------------------

elina_poly.ml elina_poly.mli elina_poly_caml.c: elina_poly.idl perlscript_caml_poly.pl 
	mkdir -p tmp
	cp $(MLAPRONDIR)/manager.idl elina_poly.idl tmp/
	cd tmp && $(CAMLIDL) -no-include -nocpp elina_poly.idl
	cp tmp/elina_poly_stubs.c elina_poly_caml.c
	$(PERL) perlscript_caml_poly.pl < tmp/elina_poly.ml >elina_poly.ml
	$(PERL) perlscript_caml_poly.pl < tmp/elina_poly.mli >elina_poly.mli

elina_oct.ml elina_oct.mli elina_oct_caml.c: elina_oct.idl perlscript_c_oct.pl perlscript_caml_oct.pl $(MLAPRONDIR)/*.idl
	mkdir -p tmp
	cp elina_oct.idl $(MLAPRONDIR)/*.idl tmp
	cd tmp && $(CAMLIDL) -no-include -nocpp -I . elina_oct.idl
	$(PERL) perlscript_c_oct.pl < tmp/elina_oct_stubs.c >elina_oct_caml.c
	$(PERL) perlscript_caml_oct.pl < tmp/elina_oct.ml >elina_oct.ml
	$(PERL) perlscript_caml_oct.pl < tmp/elina_oct.mli >elina_oct.mli
	


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

