PRGS =  baseml basemlg codeml codemlsites evolver mcmctree pamp yn00
CC = cl # cc, gcc, cl
CFLAGS = -W3 -O2
 
 
DEFINE = #
LIBS = #-lm -lM
 
 
OBJbaseml = baseml.obj eigen.obj tools.obj
OBJcodeml = codeml.obj eigen.obj tools.obj
 
all : $(PRGS)
 
baseml : $(OBJbaseml)
        $(CC) $(CFLAGS) $(OBJbaseml) $(LIBS)
basemlg : basemlg.obj tools.obj
        $(CC) $(CFLAGS) basemlg.obj tools.obj $(LIBS)
codeml : $(OBJcodeml)
        $(CC) $(CFLAGS) $(OBJcodeml) $(LIBS)
codemlsites : $(OBJcodeml)
        $(CC) -Fecodemlsites.exe -DNSSITESBATCH codeml.c $(CFLAGS) tools.obj eigen.obj $(LIBS)
evolver : evolver.obj tools.obj eigen.obj
        $(CC) $(CFLAGS) evolver.obj tools.obj eigen.obj $(LIBS)
mcmctree : mcmctree.obj tools.obj
        $(CC) $(CFLAGS) mcmctree.obj tools.obj $(LIBS)
pamp : pamp.obj tools.obj eigen.obj
        $(CC) $(CFLAGS) pamp.obj tools.obj eigen.obj $(LIBS)
yn00 : yn00.obj tools.obj
        $(CC) $(CFLAGS) yn00.obj tools.obj $(LIBS)


eigen.obj : tools.h eigen.c
        $(CC) $(CFLAGS) -c eigen.c
tools.obj : tools.h tools.c
        $(CC) $(CFLAGS) -c tools.c
baseml.obj : tools.h baseml.c treesub.c treespace.c
        $(CC) $(CFLAGS) -c baseml.c
basemlg.obj : tools.h basemlg.c treesub.c
        $(CC) $(CFLAGS) -c basemlg.c
codeml.obj : tools.h codeml.c treesub.c treespace.c
        $(CC) $(CFLAGS) -c codeml.c
evolver.obj: evolver.c treesub.c treespace.c
        $(CC) $(CFLAGS) -c evolver.c
mcmctree.obj : tools.h mcmctree.c treesub.c treespace.c
        $(CC) $(CFLAGS) -c mcmctree.c
pamp.obj : tools.h pamp.c treesub.c treespace.c
        $(CC) $(CFLAGS) -c pamp.c
yn00.obj : tools.h yn00.c 
        $(CC) $(CFLAGS) -c yn00.c
