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


include Makefile.config



all: c ocaml java

ifneq ($(HAS_OCAML),)
ocaml : 
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
	(cd apron_interface; make install)
ifneq ($(HAS_OCAML),) 
	(cd ocaml_interface; make all)
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
ifneq ($(HAS_OCAML),) 
	(cd ocaml_interface; make clean)
endif
ifneq ($(HAS_JAVA),) 
	(cd java_interface; make clean)
endif

