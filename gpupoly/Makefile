include ../Makefile.config

CC = nvcc 

CXXFLAGS = -Wcast-qual -Wswitch -Wall -Wextra -Wundef -Wcast-align -fPIC -O3 -DNDEBUG -std=c++14 -DNUM_LONGLONGRAT -DNUM_DOUBLE -fno-tree-vectorize -m64 -march=native

LIBDIR = $(ELINA_PREFIX)/lib

INCLDIR = $(ELINA_PREFIX)/include

LIBS = -lcublas

INCLUDES = -Iinclude/ -Isrc/ -I.


INSTALL = install
INSTALLd = install -d

OBJS = bindings.o maxpool2d.o affineexpr.o conv2d.o dense.o mmm.o mvm.o network.o relu.o matrix.o vector.o filters.o

SOINST = libgpupoly.so

all : libgpupoly.so

bindings.o : bindings.cu
	$(CC) -ccbin g++ -c --compiler-options "$(CXXFLAGS)" -arch=sm_61 $(INCLUDES) -o bindings.o bindings.cu $(LIBS)

maxpool2d.o : layers/maxpool2d.h layers/maxpool2d.cu
	$(CC) -ccbin g++ -c --compiler-options "$(CXXFLAGS)" -arch=sm_61 $(INCLUDES) -o maxpool2d.o layers/maxpool2d.cu $(LIBS)
	
affineexpr.o : affineexpr.h affineexpr.cu
	$(CC) -ccbin g++ -c --compiler-options "$(CXXFLAGS)" -arch=sm_61 $(INCLUDES) -o affineexpr.o affineexpr.cu $(LIBS)

conv2d.o: layers/conv2d.h layers/conv2d.cu
	$(CC) -ccbin g++ -c --compiler-options "$(CXXFLAGS)" -arch=sm_61 $(INCLUDES) -o conv2d.o layers/conv2d.cu $(LIBS)

dense.o: layers/dense.h layers/dense.cu
	$(CC) -ccbin g++ -c --compiler-options "$(CXXFLAGS)" -arch=sm_61 $(INCLUDES) -o dense.o layers/dense.cu $(LIBS)

mmm.o : mmm.cu
	$(CC) -ccbin g++ -c --compiler-options "$(CXXFLAGS)" -arch=sm_61 $(INCLUDES) -o mmm.o mmm.cu $(LIBS)

mvm.o : mvm.cu
	$(CC) -ccbin g++ -c --compiler-options "$(CXXFLAGS)" -arch=sm_61 $(INCLUDES) -o mvm.o mvm.cu $(LIBS)

network.o : network.h network.cu
	$(CC) -ccbin g++ -c --compiler-options "$(CXXFLAGS)" -arch=sm_61 $(INCLUDES) -o network.o network.cu $(LIBS)
	
relu.o: layers/relu.h layers/relu.cu
	$(CC) -ccbin g++ -c --compiler-options "$(CXXFLAGS)" -arch=sm_61 $(INCLUDES) -o relu.o layers/relu.cu $(LIBS)
	
matrix.o : matrix.h matrix.cu
	$(CC) -ccbin g++ -c --compiler-options "$(CXXFLAGS)" -arch=sm_61 $(INCLUDES) -o matrix.o matrix.cu $(LIBS)
	
vector.o : vector.h vector.cu
	$(CC) -ccbin g++ -c --compiler-options "$(CXXFLAGS)" -arch=sm_61 $(INCLUDES) -o vector.o vector.cu $(LIBS)

filters.o : filters.h filters.cu
	$(CC) -ccbin g++ -c --compiler-options "$(CXXFLAGS)" -arch=sm_61 $(INCLUDES) -o filters.o filters.cu $(LIBS)

libgpupoly.so : $(OBJS) 
	$(CC) -ccbin g++ -shared $(CC_ELINA_DYLIB) --compiler-options "$(CXXFLAGS)" -arch=sm_61 $(INCLUDES) -o $(SOINST) $(OBJS) $(LIBS)


install:
	$(INSTALLd) $(LIBDIR); \
	for i in $(SOINST); do \
		$(INSTALL) $$i $(LIBDIR); \
	done;

	

clean:
	-rm *.o
	-rm *.so





