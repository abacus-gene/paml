PRGS =  baseml basemlg codeml codemlsites evolver mcmctree pamp yn00 chi2
CC = cl # cc, gcc, cl
CFLAGS =  -O2 -Ot -Oi -Og -G6 -nologo -ML -Za -W3 -DWIN32 -DNDEBUG



DEFINE = #
LIBS = #-lm -lM

all : $(PRGS)
 
baseml : baseml.obj eigen.obj tools.obj
  $(CC) $(CFLAGS) baseml.obj eigen.obj tools.obj $(LIBS)
basemlg : basemlg.obj tools.obj
  $(CC) $(CFLAGS) basemlg.obj tools.obj $(LIBS)
codeml : codeml.obj eigen.obj tools.obj
  $(CC) $(CFLAGS) codeml.obj eigen.obj tools.obj $(LIBS)
codemlsites : codemlsites.obj eigen.obj tools.obj
  $(CC) $(CFLAGS) -Fecodemlsites.exe codemlsites.obj tools.obj eigen.obj $(LIBS)
evolver : evolver.obj tools.obj eigen.obj
  $(CC) $(CFLAGS) evolver.obj tools.obj eigen.obj $(LIBS)
mcmctree : mcmctree.obj tools.obj
  $(CC) $(CFLAGS) mcmctree.obj tools.obj $(LIBS)
pamp : pamp.obj tools.obj eigen.obj
  $(CC) $(CFLAGS) pamp.obj tools.obj eigen.obj $(LIBS)
yn00 : yn00.obj tools.obj
  $(CC) $(CFLAGS) yn00.obj tools.obj $(LIBS)
chi2 : chi2.c
  $(CC) $(CFLAGS) chi2.c $(LIBS)


eigen.obj : paml.h eigen.c
  $(CC) $(CFLAGS) -c eigen.c
tools.obj : paml.h tools.c
  $(CC) $(CFLAGS) -c tools.c
baseml.obj : paml.h baseml.c treesub.c treespace.c
  $(CC) $(CFLAGS) -c baseml.c
basemlg.obj : paml.h basemlg.c treesub.c
  $(CC) $(CFLAGS) -c basemlg.c
codeml.obj : paml.h codeml.c treesub.c treespace.c
  $(CC) $(CFLAGS) -c codeml.c
codemlsites.obj : paml.h codeml.c treesub.c treespace.c
  $(CC) $(CFLAGS) -c -Focodemlsites.obj -DNSSITESBATCH codeml.c
evolver.obj: evolver.c treesub.c treespace.c
  $(CC) $(CFLAGS) -c evolver.c
mcmctree.obj : paml.h mcmctree.c treesub.c treespace.c
  $(CC) $(CFLAGS) -c mcmctree.c
pamp.obj : paml.h pamp.c treesub.c treespace.c
  $(CC) $(CFLAGS) -c pamp.c
yn00.obj : paml.h yn00.c 
  $(CC) $(CFLAGS) -c yn00.c
