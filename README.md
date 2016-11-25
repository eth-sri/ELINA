# ELINA
ETH Library for Numerical Analyis (ELINA) contains optimized implementations of popular numerical abstract domains such as Polyhedra and Octagons for static analysis. It builds on top of APRON which is a popular library for static analysis. 

The library uses improved algorithms, online decomposition as well as state of the art performance optimizations from linear algebra such as vectorization, locality of reference, scalar replacement etc. to significantly improve the performance of static analysis with the numerical domains.

#Requirements:
  The installation is preferable for Linux 64-bit. Make sure you have latest version of gcc and g++.

  Install the following libraries.

    1. GNU m4:
		sudo apt-get install m4

    2. GMP:
		a. Download the tar file from https://gmplib.org/.
		b. Extract the source.
		c. Go to the gmp folder and run:
			./configure --enable-cxx
			make 
			make check
			sudo make install
		d. This will install the gmp library in "/usr/local" folder.

    3. GNU MPFR:
		a. Download the tar file from http://www.mpfr.org/
		b. Extract the source.
		c. Go to the mpfr folder and run:
			./configure
			make
			make check
			sudo make install
		d. This will install the mpfr library in "/usr/local" folder.

    4. APRON:
		a. Download source from http://apron.cri.ensmp.fr/library/
        	b. Go to the "apron" folder.
		c. Install the library as per README file. 
		d. Specify "/usr/local" for gmp and mpfr library paths in "Makefile.config".
  
#Compiling and installing:
    Copy the "partitions_api", "elina_oct", "elina_poly" folder into the "apron" directory. The "partitions_api" needs to be compiled and installed first by running "make" and then "sudo make install" in the "partitions_api" folder.
    To compile and install the Octagon domain do the following:
		1. Check if your machine supports SSE or AVX as follows:
			a. Run "cat /proc/cpuinfo | grep "sse\|avx"".
			b. Check for strings "sse", "avx" in output.
		2. To compile the library: 
			a. If your computer supports AVX 
			   run "make IS_VECTOR=-DVECTOR" to compile the source. 
			   This will use AVX vectorized operators.
			b. Else if your computer supports SSE 
			   run "make IS_VECTOR=-DVECTOR IS_SSE=-DSSE" to compile the source. 
			   This will use SSE vectorized operators.
			c. Else run "make". This will use scalar operators.
		3. For 32-bit systems, disable the required architecture flags (m64, march=native etc. ) and DTIMING in DFLAGS inside the Makefile
		4. The Octagon domain can be installed by running "sudo make install".
      
     The Polyhedra domain does not require vector intrinsics and can be compiled by running "make" in the "elina_poly" folder. For 32-bit systems, make sure to disable the architecture flags in the Makefile. The generated shared object can be installed by running "sudo make install".

    
#How to Use in Static Analyzer
  The library can be used directly from C. Besides this, we also provide interfaces for C++ and Java also.
  These interfaces can be installed as follows:

  Java:
	
      Copy the source files in "java_interface" directory into "japron/apron" directory
      Replace the Makefile in "japron" directory with the Makefile in "java_interface" directory.
      run "make" in "japron" directory.
      run "sudo make install" to install the updated "libjapron.so" file.
      Initialize the APRON Manager for Octagon as:
        man = new OptOctagon();
      and for Polyhedra as:
	man = new OptPoly(false);
      
  C++:
      
      Copy the source files in "C++ interface" directory into "apronxx" directory
      Replace the Makefile in "apronxx" directory with the Makefile in "C++ interface" directory.
      run "make" in "apronxx" directory.
      run "sudo make install" to install the updated "libapronxx.so" file.
      Initialize the APRON Manager for Octagon as:
        man = new opt_oct_manager();
      and for Polyhedra as:
	man = new opt_pk_manager();
