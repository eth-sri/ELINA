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


include Makefile.config



all: c ocaml java

ifneq ($(HAS_OCAML),)
ocaml : c
	(cd ocaml_interface; make all)
else
ocaml : 
endif

ifneq ($(HAS_JAVA),)
java : 
	(cd java_interface; make all)
else
java : 
endif

c:
ifeq ($(IS_APRON),)
	(cd elina_auxiliary; make all)
endif
	(cd elina_linearize; make all)
	(cd partitions_api; make all)
	(cd elina_oct; make all)
	(cd elina_poly; make all)
	(cd elina_zones; make all)
	(cd elina_zonotope; make all) 
	(cd zonoml; make all)
	(cd fppoly; make all)


install:
ifeq ($(IS_APRON),)
	(cd elina_auxiliary; make install)
endif
	(cd elina_linearize; make install)
	(cd partitions_api; make install)
	(cd elina_oct; make install)
	(cd elina_poly; make install)
	(cd elina_zones; make install)
	(cd elina_zonotope; make install)
	(cd zonoml; make install)
	(cd fppoly; make install)
	(cd apron_interface; make install)
ifneq ($(HAS_OCAML),) 
	(cd ocaml_interface; make install)
ifneq ($(OCAMLFIND),)
	$(OCAMLFIND) remove elina
	$(OCAMLFIND) install elina ocaml_interface/META ocaml_interface/dllelina_poly_caml.so ocaml_interface/elina_poly.a ocaml_interface/elina_poly.cma ocaml_interface/elina_poly.cmi ocaml_interface/elina_poly.cmo ocaml_interface/elina_poly.cmx ocaml_interface/elina_poly.cmxa ocaml_interface/elina_poly.idl ocaml_interface/elina_poly.ml ocaml_interface/elina_poly.mli ocaml_interface/elina_poly.o ocaml_interface/elina_poly_caml.c ocaml_interface/elina_poly_caml.o ocaml_interface/libelina_poly_caml.a
endif
endif
ifneq ($(HAS_JAVA),) 
	(cd java_interface; make all)
endif

clean:
ifeq ($(IS_APRON),)
	(cd elina_auxiliary; make clean)
endif
	(cd elina_linearize; make clean)
	(cd partitions_api; make clean)
	(cd elina_oct; make clean)
	(cd elina_poly; make clean)
	(cd elina_zones; make clean)
	(cd elina_zonotope; make clean)
	(cd zonoml; make clean)
	(cd fppoly; make clean)
ifneq ($(HAS_OCAML),) 
	(cd ocaml_interface; make clean)
endif
ifneq ($(HAS_JAVA),) 
	(cd java_interface; make clean)
endif

