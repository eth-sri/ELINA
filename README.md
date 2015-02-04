# Octagon
An Optimized implementation of the Octagon Abstract Domain

#Requirements:
  Make sure you have latest version of APRON library installed. The APRON library can be downloaded from
    http://apron.cri.ensmp.fr/library/
  
#Compiling:
  Download the source and follow the instructions below.
    Copy the "src" folder into the directory containing APRON so that the path to new folder is: APRON_PATH/src
    Go to the "APRON_PATH/src" folder and open the Makefile
    If your machine supports Intel's AVX intrinsics, and you would like to use vectorized operators for dense type,
      specify -DVECTOR for IS_VECTOR and "-m64", "-march=native" for AFLAGS in Makefile
      otherwise comment it out
      specify -DTHRESHOLD=value, for example -DTHRESHOLD=0.75 in DFLAGS. The Threshold lies between 0 and 1 and controls       switching between dense and decomposed types (operators). The analysis will switch from decomposed to dense as soon       as sparsity is below the threshold. A larger value of threshold favors dense whereas smaller favors decomposed           types(operators).
      Run "make", it will generate "liboptoct.so" and "liblinkedlistapi.so" libraries.
      
  #Installing:
    Run "sudo make install" to install the libraries into the directory specified in APRON's "Makefile.config" file.

