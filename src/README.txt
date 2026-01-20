Notes for compiling PAML on UNIX systems, including MAC OS X

Ziheng Yang (z.yang@ucl.ac.uk)
Last updated, 20 January 2026


Copyright notice and disclaimer

The software package is provided "as is" without warranty of any
kind. In no event shall the author or his employer be held responsible
for any damage resulting from the use of this software, including but
not limited to the frustration that you may experience in using the
package.  The program package, including source codes, example data
sets, executables, and this documentation, is maintained by Ziheng
Yang and distributed under the GNU GPL v3.


Instructions for compiling 
==========================

Method I

Use the Makefile.  The default Makefile is for UNIX/Linux/OSX.  

   make

   cp baseml basemlg codeml evolver pamp yn00 mcmctree chi2 ..
   rm *.o
   cd ..

On some systems you might have to edit the Makefile and change a few
flags at the beginning of the file.

To compile using MSVC++, use Makefile.MSVC:
   nmake -f Makefile.MSVC


Method II

You can also compile the programs from the command line.  Here are the
commands for the cc and gcc compilers.  You might have to recombine
the different choices (cc vs. gcc, -fast vs. -O2 or -O3, and with or
without -lm).  Make sure you turn on some optimization options (-O2,
-O3, -fast, etc.) as otherwise the code can be several times slower.
Below are a few possibilities.

(2a) MAC OS X Developer's Toolkit

cc -O2 -o baseml baseml.c tools.c -lm
cc -O2 -o basemlg basemlg.c tools.c -lm
cc -O2 -o codeml codeml.c tools.c -lm
cc -O2 -o pamp pamp.c tools.c -lm
cc -O2 -o mcmctree mcmctree.c tools.c -lm
cc -O2 -o infinitesites -D INFINITESITES mcmctree.c tools.c -lm
cc -O2 -o evolver evolver.c tools.c -lm
cc -O2 -o yn00 yn00.c tools.c -lm
cc -O2 -o chi2 chi2.c -lm
cc -O2 -o ds ds.c tools.c -lm


(2b) gcc compiler

gcc -O3 -o baseml baseml.c tools.c
gcc -O3 -o basemlg basemlg.c tools.c
gcc -O3 -o codeml codeml.c tools.c
gcc -O3 -o pamp pamp.c tools.c
gcc -O3 -o mcmctree mcmctree.c tools.c
gcc -O3 -o infinitesites -D INFINITESITES mcmctree.c tools.c -lm
gcc -O3 -o evolver evolver.c tools.c
gcc -O3 -o yn00 yn00.c tools.c 
gcc -O3 -o chi2 chi2.c 
gcc -o ds -O3 ds.c tools.c -lm

> NOTE: These commands should work on the Windows Command Prompt
>       You could write the extension ".exe" after the PAML 
>       program but, oftentimes, gcc adds this automatically.

(2c) UNIX cc compiler

cc -fast -o baseml baseml.c tools.c -lm
cc -fast -o basemlg basemlg.c tools.c -lm
cc -fast -o codeml codeml.c tools.c -lm
cc -fast -o pamp pamp.c tools.c -lm
cc -fast -o mcmctree mcmctree.c tools.c -lm
cc -fast -o infinitesites -D INFINITESITES mcmctree.c tools.c -lm
cc -fast -o evolver evolver.c tools.c -lm
cc -fast -o yn00 yn00.c tools.c -lm
cc -fast -o chi2 chi2.c -lm
cc -fast -o ds ds.c tools.c -lm


If you have the root or administrator, you can copy the executables
into the folder /usr/local/bin/, so that the programs are available to
everyone with an account on the system.

scp  baseml basemlg codeml pamp mcmctree evolver yn00 chi2 infinitesites ds  /usr/local/bin/


// End of file
