PRGS =  baseml codeml basemlg mcmctree pamp evolver yn00 chi2 
CC = cc # cc, gcc, cl

CFLAGS = -O3
#CFLAGS = -fast

LIBS = -lm # -lM

all : $(PRGS)

baseml : baseml.c tools.c treesub.c treespace.c paml.h 
	$(CC) $(CFLAGS) -o $@ baseml.c tools.c $(LIBS)
basemlg : basemlg.c  tools.c treesub.c treespace.c paml.h 
	$(CC) $(CFLAGS) -o $@ basemlg.c tools.c $(LIBS)
codeml : codeml.c  tools.c treesub.c treespace.c paml.h 
	$(CC) $(CFLAGS) -o $@ codeml.c tools.c $(LIBS)
evolver : evolver.c  tools.c treesub.c treespace.c paml.h 
	$(CC) $(CFLAGS) -o $@ evolver.c tools.c $(LIBS)
pamp : pamp.c  tools.c treesub.c treespace.c paml.h 
	$(CC) $(CFLAGS) -o $@ pamp.c tools.c $(LIBS)
mcmctree : mcmctree.c  tools.c treesub.c treespace.c paml.h 
	$(CC) $(CFLAGS) -o $@ mcmctree.c tools.c $(LIBS)
	$(CC) $(CFLAGS) -o infinitesites -D INFINITESITES mcmctree.c tools.c $(LIBS)
yn00: yn00.c tools.c paml.h
	$(CC) $(CFLAGS) -o $@ yn00.c tools.c $(LIBS)
chi2 : chi2.c
	$(CC) $(CFLAGS) -o $@ chi2.c $(LIBS)
