PRGS =  baseml codeml basemlg pamp evolver mcmctree # mljc3s dsdn
CC = cc # cc, gcc, cl

CFLAGS = -w0 -warnprotos -newc -std -fast

#CFLAGS = -w0 -warnprotos -newc -std -fast
#CFLAGS = -g -w0 -warnprotos -newc -std -edit1 #
#CFLAGS = -g -w0 -warnprotos -newc -std  #


DEFINE = #
LIBS = -lm # -lM

OBJbaseml = baseml.o eigen.o tools.o
OBJcodeml = codeml.o eigen.o tools.o

all : $(PRGS)

dsdn: dsdn.c tools.o eigen.o
	$(CC) $(CFLAGS) -o $@ dsdn.c tools.o eigen.o $(LIBS)

mljc3s: mljc3s.c tools.o
	$(CC) $(CFLAGS) -o $@ mljc3s.c tools.o $(LIBS)

baseml : $(OBJbaseml)
	$(CC) $(CFLAGS) -o $@ $(OBJbaseml) $(LIBS)
basemlg : basemlg.o tools.o
	$(CC) $(CFLAGS) -o $@ basemlg.o tools.o $(LIBS)
codeml : $(OBJcodeml)
	$(CC) $(CFLAGS) -o $@ $(OBJcodeml) $(LIBS)
evolver : evolver.o tools.o 
	$(CC) $(CFLAGS) -o $@ tools.o evolver.o $(LIBS)
pamp : pamp.o tools.o eigen.o
	$(CC) $(CFLAGS) -o $@ pamp.o tools.o eigen.o $(LIBS)
mcmctree : mcmctree.o tools.o
	$(CC) $(CFLAGS) -o $@ mcmctree.o tools.o $(LIBS)

eigen.o : tools.h eigen.c
	$(CC) $(CFLAGS) -c eigen.c
tools.o : tools.h tools.c
	$(CC) $(CFLAGS) -c tools.c
baseml.o : tools.h baseml.c treesub.c treespace.c
	$(CC) $(CFLAGS) -c baseml.c
basemlg.o : tools.h basemlg.c treesub.c
	$(CC) $(CFLAGS) -c basemlg.c
codeml.o : tools.h codeml.c treesub.c treespace.c
	$(CC) $(CFLAGS) -c codeml.c
evolver.o: evolver.c treesub.c treespace.c
	$(CC) $(CFLAGS) -c evolver.c
pamp.o : tools.h pamp.c treesub.c treespace.c
	$(CC) $(CFLAGS) -c pamp.c
mcmctree.o : tools.h mcmctree.c treesub.c treespace.c
	$(CC) $(CFLAGS) -c mcmctree.c
