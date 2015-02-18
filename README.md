# Octagon
An Optimized Library for the Octagon Abstract Domain

#Requirements:
  Make sure you have latest version of APRON library installed. The APRON library can be downloaded from
    http://apron.cri.ensmp.fr/library/
  
#Compiling:

    Download the source and follow the instructions below.

    Copy the "optoctagons" folder into the directory containing APRON so that the path to new folder is: APRON_PATH/optoctagons.

    Go to the "APRON_PATH/optoctagons" folder and open the Makefile.
   
   Go to the "apron/optoctagons" directory in terminal
	
		1. Check if your machine supports SSE or AVX as follows:
			a. Run "cat /proc/cpuinfo | grep "sse\|avx"".
			b. Check for strings "sse", "avx" in output.

		2. To compile the library: 
			a. If your computer supports AVX run "make IS_VECTOR=-DVECTOR" to compile the source. This will use AVX vectorized operators for dense type.
			b. Else if your computer supports SSE run "make IS_VECTOR=-DVECTOR IS_SSE=-DSSE" to compile the source. This will use SSE vectorized operators for dense type.
			c. Else run "make". This will use scalar operators for dense type.
      
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
