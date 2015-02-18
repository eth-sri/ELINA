# OptOctagons
OptOctagons is an optimized library for static program analysis with the Octagon numerical domain. It builds on top of APRON which is a popular library for static analysis. 

The library uses improved algorithms, online decomposition of octagons as well as state of the art performance optimizations from linear algebra such as vectorization, locality of reference, scalar replacement etc. to significantly improve the performance of static analysis with the Octagon domain.

#Requirements:
  Make sure you have latest version of APRON library installed. The APRON library can be downloaded from
    http://apron.cri.ensmp.fr/library/
  
#Compiling:
    Copy the "optoctagons" folder into the APRON directory.
    Go to the "APRON_PATH/optoctagons" folder in terminal.
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
      
#Installing:
    Specify the install directory in APRON's "Makefile.config" file.
    Run "sudo make install"
    
#How to Use in Static Analyzer

  Java:
	
      Copy the files in "java_interface" directory into APRON_PATH/japron/apron
      Replace the Makefile in APRON_PATH/japron directory with the Makefile in "java_interface" directory.
      run "make" in the APRON_PATH/japron
      run "sudo make install" to install the updated "libjapron.so" file.
      Initialize the APRON Manager as:
        man = new OptOctagon();
      
  C++:
      
      Copy the files in "C++ interface" directory into APRON_PATH/apronxx
      Replace the Makefile in APRON_PATH/apronxx directory with the Makefile in "C++ interface" directory.
      run "make" in the APRON_PATH/japron
      run "sudo make install" to install the updated "libjapron.so" file.
      Initialize the APRON Manager as:
        man = new opt_oct_manager();
