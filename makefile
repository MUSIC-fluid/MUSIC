/SHELL=/bin/sh
# This file contains a set of rules used by the "make"
#   command.  This makefile tells "make" how the
#   executable program $(COMMAND) should be created
#   from the source files $(SRCS) via the object
#   files $(OBJS).
 
# This file also tells make what to do if arguments
#   follow the "make" command.
#
# To remove the OBJS files; issue the command:
#        "make clean"
#
# To create a tar file with name $(COMMAND).tar
#  containing this makefile and the SRCS files:
#  issue the command:
#        "make tar"
 
SRCS= \
main.cpp eos.cpp evolve.cpp grid.cpp  \
init.cpp reconst.cpp freeze.cpp minmod.cpp\
random.cpp glauber.cpp util.cpp advance.cpp u_derivative.cpp dissipative.cpp
  
OBJS= \
grid.o  eos.o    evolve.o     \
init.o    reconst.o  freeze.o minmod.o  \
glauber.o advance.o u_derivative.o \
dissipative.o random.o\
util.o main.o 

 
CFT= mpiCC
FFLAGS= -O3 -Wno-write-strings
O3FFLAGS=
LDFLAGS=
LIBS= -lm -lgsl -lblas
COMMAND=  mpihydro
 
 
$(COMMAND): $(OBJS)
	$(CFT) $(FFLAGS) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(LIBS)

minmod.o : minmod.cpp minmod.h 
	$(CFT) $(FFLAGS)  -c minmod.cpp -o minmod.o

advance.o : advance.cpp advance.h util.h dissipative.h minmod.h grid.h
	$(CFT) $(FFLAGS)  -c advance.cpp -o advance.o

dissipative.o : dissipative.cpp dissipative.h minmod.h util.h grid.h
	$(CFT) $(FFLAGS)  -c dissipative.cpp -o dissipative.o

u_derivative.o : u_derivative.cpp u_derivative.h util.h
	$(CFT) $(FFLAGS)  -c u_derivative.cpp -o u_derivative.o

eos.o : eos.cpp eos.h util.h 
	$(CFT) $(FFLAGS)  -c eos.cpp -o eos.o

random.o : random.cpp random.h util.h 
	$(CFT) $(FFLAGS)  -c random.cpp -o random.o
 
evolve.o : evolve.cpp data.h eos.h evolve.h grid.h reconst.h advance.h util.h grid.h
	$(CFT) $(FFLAGS)  -c evolve.cpp -o evolve.o
 
grid.o : grid.cpp grid.h util.h 
	$(CFT) $(FFLAGS)  -c grid.cpp -o grid.o
 
init.o : init.cpp eos.h grid.h init.h util.h 
	$(CFT) $(FFLAGS)  -c init.cpp -o init.o

int.o : int.cpp util.h 
	$(CFT) $(FFLAGS)  -c int.cpp -o int.o
 
reconst.o : reconst.cpp data.h eos.h grid.h reconst.h util.h 
	$(CFT) $(FFLAGS)  -c reconst.cpp -o reconst.o
 
util.o : util.cpp util.h 
	$(CFT) $(FFLAGS)  -c util.cpp -o util.o

glauber.o : glauber.cpp util.h glauber.h
	$(CFT) $(FFLAGS)  -c glauber.cpp -o glauber.o

freeze.o : freeze.cpp data.h eos.h grid.h util.h freeze.h int.h
	$(CFT) $(FFLAGS)  -c freeze.cpp -o freeze.o

main.o : main.cpp data.h eos.h evolve.h grid.h init.h util.h 
	$(CFT) $(FFLAGS)  -c main.cpp -o main.o

clean:
	rm -f $(OBJS)
 
tar:
	tar cf - makefile eos.cpp evolve.cpp grid.cpp init.cpp main.cpp reconst.cpp util.cpp glauber.cpp eos.h evolve.h grid.h init.h util.h glauber.h data.h freeze.cpp freeze.h input > $(COMMAND).tar

compress:	tar
	compress $(COMMAND).tar
