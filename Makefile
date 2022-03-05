## CPLEX
SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio1271/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio1271/concert
#CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio_Preview124/cplex
#CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio_Preview124/concert
CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD
COPT  = -m64 -fPIC -fno-strict-aliasing 
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -lm -pthread
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
EXDIR         = $(CPLEXDIR)/examples
EXINC         = $(EXDIR)/include
EXDATA        = $(EXDIR)/data
EXSRCC        = $(EXDIR)/src/c
EXSRCCX       = $(EXDIR)/src/c_x
## /CPLEX

CC       = gcc  -DPUC -O3
#CFLAGS  = -x03
CFLAGS  = -g $(COPT)  -I$(CPLEXINCDIR)
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 


EXECS   = pm_path 

all: $(EXECS)

clean:
	/bin/rm -f *.o;

bibrand.o: bibrand.c
	$(CC)  $(CFLAGS) -c bibrand.c -o bibrand.o

2pndp.o: 2pndp.c 
	$(CC)  $(CFLAGS) -c 2pndp.c -o 2pndp.o 

best.o: best.c
	$(CC)  $(CFLAGS) -c best.c -o best.o

mine.o: mine.c
	$(CC)  $(CFLAGS) -c mine.c -o mine.o

grasp-sem.o: grasp-sem.c
	$(CC)  $(CFLAGS) -c grasp-sem.c -o grasp-sem.o

grasp-sem: grasp-sem.o 2pndp.o bibrand.o best.o mine.o
	$(CC) $(CFLAGS) -o grasp-sem grasp-sem.o 2pndp.o bibrand.o best.o mine.o -lm

grasp-prbf-alvo.o: grasp-prbf-alvo.c
	$(CC)  $(CFLAGS) -c grasp-prbf-alvo.c -o grasp-prbf-alvo.o

grasp-prbf-alvo: grasp-prbf-alvo.o 2pndp.o bibrand.o best.o mine.o
	$(CC) $(CFLAGS) -o grasp-prbf-alvo grasp-prbf-alvo.o 2pndp.o bibrand.o best.o mine.o -lm

grasp-prbf-sem.o: grasp-prbf-sem.c
	$(CC)  $(CFLAGS) -c grasp-prbf-sem.c -o grasp-prbf-sem.o

grasp-prbf-sem: grasp-prbf-sem.o 2pndp.o bibrand.o best.o mine.o
	$(CC) $(CFLAGS) -o grasp-prbf-sem grasp-prbf-sem.o 2pndp.o bibrand.o best.o mine.o -lm

grasp-prbf-time.o: gprbf-time.c
	$(CC)  $(CFLAGS) -c gprbf-time.c -o grasp-prbf-time.o

grasp-prbf-time: grasp-prbf-time.o 2pndp.o bibrand.o best.o mine.o
	$(CC) $(CFLAGS) -o grasp-prbf-time grasp-prbf-time.o 2pndp.o bibrand.o best.o mine.o -lm

grasp.o: grasp.c
	$(CC)  $(CFLAGS) -c grasp.c -o grasp.o

grasp: grasp.o 2pndp.o bibrand.o best.o mine.o
	$(CC) $(CFLAGS) -o grasp grasp.o 2pndp.o bibrand.o best.o mine.o -lm

pm_2star.o: pm_2star.c
	$(CC) -c $(CFLAGS) pm_2star.c -o pm_2star.o

pm_2star: pm_2star.o 2pndp.o bibrand.o best.o mine.o
	$(CC) $(CFLAGS) -o pm_2star pm_2star.o 2pndp.o bibrand.o best.o mine.o $(CLNFLAGS)

pm_path.o: pm_path.c
	$(CC) -c $(CFLAGS) pm_path.c -o pm_path.o

pm_path: pm_path.o 2pndp.o bibrand.o best.o mine.o
	$(CC) $(CFLAGS) -o pm_path pm_path.o 2pndp.o bibrand.o best.o mine.o $(CLNFLAGS)

