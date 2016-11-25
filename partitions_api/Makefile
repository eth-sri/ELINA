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

DFLAGS = -g

AFLAGS := -D_GNU_SOURCE -pthread -m64 -march=native -fno-tree-vectorize 

SOINST = libpartitions.so

OBJS = comp_list.o array_comp_list.o intersection.o union.o extract.o 

PREFIX = $(APRON_PREFIX)

LIBDIR = $(PREFIX)/lib

INCLDIR = $(PREFIX)/include

PARTITIONSH = comp_list.h

all: libpartitions.so

comp_list.o: comp_list.h comp_list.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(AFLAGS) -o comp_list.o comp_list.c $(LIBS)

array_comp_list.o: comp_list.h array_comp_list.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(AFLAGS) -o array_comp_list.o array_comp_list.c $(LIBS)

intersection.o: comp_list.h intersection.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(AFLAGS) -o intersection.o intersection.c $(LIBS)

union.o: comp_list.h union.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(AFLAGS) -o union.o union.c $(LIBS)


#closure_binary.o: closure_binary.h closure_binary.c
#	$(CC) -c $(CFLAGS) $(DFLAGS) $(AFLAGS) -o closure_binary.o closure_binary.c $(LIBS)

extract.o: comp_list.h extract.c
	$(CC) -c $(CFLAGS) $(DFLAGS) $(AFLAGS) -o extract.o extract.c $(LIBS)

libpartitions.so: $(OBJS)
	$(CC) -shared $(CFLAGS) $(DFLAGS) -o $(SOINST) $(OBJS) $(LIBS)

install: 
	$(INSTALLd) $(LIBDIR); \
	$(INSTALL) $(SOINST) $(LIBDIR);
	$(INSTALLd) $(INCLDIR); \
	for i in $(PARTITIONSH); do \
		$(INSTALL) $$i $(INCLDIR); \
	done; 


clean:
	rm $(OBJS) libpartitions.so
