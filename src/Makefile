PRGS =  baseml codeml basemlg mcmctree pamp evolver yn00 chi2
CC = gcc
CFLAGS = -O3 -Wall -Wno-unused-variable -Wno-unused-result

LIBS = -lm

all : $(PRGS)

baseml : baseml.o tools.o treesub.c treespace.c paml.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ baseml.o tools.o $(LIBS)
basemlg: basemlg.o tools.o treesub.c treespace.c paml.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ basemlg.o tools.o $(LIBS)
codeml: codeml.o tools.o treesub.c treespace.c paml.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ codeml.o tools.o $(LIBS)
evolver: evolver.o tools.o treesub.c treespace.c paml.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ evolver.o tools.o $(LIBS)
pamp: pamp.o tools.o treesub.c treespace.c paml.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ pamp.o tools.o $(LIBS)
mcmctree: mcmctree.o tools.o treesub.c treespace.c paml.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ mcmctree.c tools.o $(LIBS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o infinitesites -D INFINITESITES mcmctree.c tools.o $(LIBS)
yn00: yn00.o tools.o paml.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ yn00.o tools.o $(LIBS)
chi2: chi2.o
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ chi2.c $(LIBS)

tools.o : paml.h tools.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -c tools.c
baseml.o : paml.h baseml.c treesub.c treespace.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -c baseml.c
basemlg.o : paml.h basemlg.c treesub.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -c basemlg.c
codeml.o : paml.h codeml.c treesub.c treespace.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -c codeml.c
evolver.o : evolver.c treesub.c treespace.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -c evolver.c
mcmctree.o : paml.h mcmctree.c treesub.c treespace.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -c mcmctree.c
pamp.o : paml.h pamp.c treesub.c treespace.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -c pamp.c
yn00.o : paml.h yn00.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -c yn00.c

clean :
	-rm *.o $(PRGS)
