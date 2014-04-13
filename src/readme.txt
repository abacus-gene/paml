Notes for compiling PAML on UNIX

Ziheng Yang (z.yang@ucl.ac.uk)
15 July 2001


(1) A Makefile is provided.  In most UNIX systems, you can simply type the command 

   make

and it will work.  In some systems you might need to change a few
flags at the beginning of the Makefile.  If you like, you can then move the
executables one level up into the paml directory and out of the src.


(2) You can also compile the programs from the command line.  Here
are the commands for the cc and gcc compiler.  You might have to
recombine the different choices (cc vs. gcc, -fast vs. -O2 or -O3, and
with or without -lm).

(a)

cc -c -fast tools.c
cc -c -fast eigen.c
cc -o baseml -fast baseml.c tools.o eigen.o -lm
cc -o codeml -fast  codeml.c tools.o eigen.o -lm
cc -o pamp -fast pamp.c tools.o eigen.o -lm
cc -o mcmctree -fast mcmctree.c tools.o -lm
cc -o evolver -fast evolver.c eigen.o tools.o -lm
cc -o yn00 -fast yn00.c tools.o -lm


(b) 

gcc -c -O2 tools.c
gcc -c -O2 eigen.c
gcc -o baseml -O2 baseml.c tools.o eigen.o 
gcc -o basemlg -O2 basemlg.c tools.o
gcc -o codeml -O2  codeml.c tools.o eigen.o
gcc -o pamp -O2 pamp.c tools.o eigen.o
gcc -o mcmctree -O2 mcmctree.c tools.o 
gcc -o evolver -O2 evolver.c tools.o eigen.o 
gcc -o yn00 -O2 yn00.c tools.o 


(c)

gcc -c -O2 tools.c eigen.c
gcc -o baseml -O2 baseml.c tools.o eigen.o -lm
gcc -o basemlg -O2 basemlg.c tools.o -lm
gcc -o codeml -O2  codeml.c tools.o eigen.o -lm
gcc -o pamp -O2 pamp.c tools.o eigen.o -lm
gcc -o mcmctree -O2 mcmctree.c tools.o -lm
gcc -o evolver -O2 evolver.c tools.o eigen.o -lm
gcc -o yn00 -O2 yn00.c tools.o -lm

// end
