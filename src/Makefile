PRGS =  baseml codeml basemlg mcmctree pamp evolver yn00 chi2
CC = cc # cc, gcc, cl
CFLAGS = -O3 -Wall -Wno-unused-result
#CC = icc
#CFLAGS = -fast -Wall

LIBS = -lm # -lM

all : $(PRGS)

baseml : baseml.o tools.o treesub.c treespace.c paml.h
	$(CC) $(CFLAGS) -o $@ baseml.o tools.o $(LIBS)
basemlg : basemlg.o tools.o treesub.c treespace.c paml.h
	$(CC) $(CFLAGS) -o $@ basemlg.o tools.o $(LIBS)
codeml : codeml.o tools.o treesub.c treespace.c paml.h
	$(CC) $(CFLAGS) -o $@ codeml.o tools.o $(LIBS)
evolver : evolver.o tools.o treesub.c treespace.c paml.h
	$(CC) $(CFLAGS) -o $@ evolver.o tools.o $(LIBS)
pamp : pamp.o tools.o treesub.c treespace.c paml.h
	$(CC) $(CFLAGS) -o $@ pamp.o tools.o $(LIBS)
mcmctree : mcmctree.o tools.o treesub.c treespace.c paml.h
	$(CC) $(CFLAGS) -o $@ mcmctree.c tools.o $(LIBS)
	$(CC) $(CFLAGS) -o infinitesites -D INFINITESITES mcmctree.c tools.o $(LIBS)
yn00: yn00.o tools.o paml.h
	$(CC) $(CFLAGS) -o $@ yn00.o tools.o $(LIBS)
chi2 : chi2.o
	$(CC) $(CFLAGS) -o $@ chi2.c $(LIBS)

tools.o : paml.h tools.c
	$(CC) $(CFLAGS) -c tools.c
baseml.o : paml.h baseml.c treesub.c treespace.c
	$(CC) $(CFLAGS) -c baseml.c
basemlg.o : paml.h basemlg.c treesub.c
	$(CC) $(CFLAGS) -c basemlg.c
codeml.o : paml.h codeml.c treesub.c treespace.c
	$(CC) $(CFLAGS) -c codeml.c
evolver.o: evolver.c treesub.c treespace.c
	$(CC) $(CFLAGS) -c evolver.c
mcmctree.o : paml.h mcmctree.c treesub.c treespace.c
	$(CC) $(CFLAGS) -c mcmctree.c
pamp.o : paml.h pamp.c treesub.c treespace.c
	$(CC) $(CFLAGS) -c pamp.c
yn00.o : paml.h yn00.c
	$(CC) $(CFLAGS) -c yn00.c

clean :
	-rm *.o $(PRGS)
