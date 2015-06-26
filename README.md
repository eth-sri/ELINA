# OptOctagons
OptOctagons is an optimized library for static program analysis with the Octagon numerical domain. It builds on top of APRON which is a popular library for static analysis. 

The library uses improved algorithms, online decomposition of octagons as well as state of the art performance optimizations from linear algebra such as vectorization, locality of reference, scalar replacement etc. to significantly improve the performance of static analysis with the Octagon domain.

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
  
#Compiling:
    Copy the "optoctagons" folder into the "apron" directory.
    Go to the new "optoctagons" folder in terminal.
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
		3. For 32-bit systems, disable the required architecture flags (m64, march=native etc) and DTIMING in DFLAGS inside the Makefile
      
#Installing:
    Specify the install directory in APRON's "Makefile.config" file.
    Run "sudo make install"
    
#How to Use in Static Analyzer
  The library can be used directly from C. Besides this, we also provide interfaces for C++ and Java also.
  These interfaces can be installed as follows:

  Java:
	
      Copy the source files in "java_interface" directory into "japron/apron" directory
      Replace the Makefile in "japron" directory with the Makefile in "java_interface" directory.
      run "make" in "japron" directory.
      run "sudo make install" to install the updated "libjapron.so" file.
      Initialize the APRON Manager as:
        man = new OptOctagon();
      
  C++:
      
      Copy the source files in "C++ interface" directory into "apronxx" directory
      Replace the Makefile in "apronxx" directory with the Makefile in "C++ interface" directory.
      run "make" in "apronxx" directory.
      run "sudo make install" to install the updated "libapronxx.so" file.
      Initialize the APRON Manager as:
        man = new opt_oct_manager();
