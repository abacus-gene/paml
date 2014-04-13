cc -c -fast tools.c
cc -c -fast eigen.c
cc -o baseml -fast baseml.c tools.o eigen.o -lm
cc -o codeml -fast  codeml.c tools.o eigen.o -lm
cc -o pamp -fast pamp.c tools.o eigen.o -lm
cc -o mcmctree -fast mcmctree.c tools.o -lm
cc -o evolver -fast evolver.c tools.o -lm
cc -o rell -fast rell.c tools.o -lm
