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
# ***Looks like this is broken***
#
# A better way:
# If you have a clone of the git repository, you
# can use 'git archive' with a command like:
# git archive --format=tar --prefix=MUSIC/ HEAD | gzip > MUSIC.tgz
 
# SRCS= \
# main.cpp eos.cpp evolve.cpp grid.cpp  \
# init.cpp reconst.cpp freeze.cpp freeze_pseudo.cpp minmod.cpp\
# random.cpp glauber.cpp util.cpp advance.cpp u_derivative.cpp dissipative.cpp
  
OBJS= \
grid.o  eos.o    evolve.o     \
init.o    reconst.o  freeze.o minmod.o  \
glauber.o advance.o u_derivative.o \
dissipative.o random.o\
util.o main.o freeze_pseudo.o

 
CXX= mpic++
CXXFLAGS=  -O3
LIBS= -L/software/libraries/GSL/1.15/lib -lm -lgsl -lgslcblas
COMMAND=  mpihydro
 
 
$(COMMAND): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(LIBS)

minmod.o : minmod.cpp minmod.h data.h

advance.o : advance.cpp advance.h util.h dissipative.h minmod.h grid.h data.h eos.h reconst.h u_derivative.h

dissipative.o : dissipative.cpp dissipative.h minmod.h util.h grid.h data.h eos.h

u_derivative.o : u_derivative.cpp u_derivative.h util.h data.h grid.h eos.h

eos.o : eos.cpp eos.h util.h data.h

random.o : random.cpp random.h util.h 

evolve.o : evolve.cpp data.h eos.h evolve.h grid.h reconst.h advance.h util.h dissipative.h minmod.h u_derivative.h

grid.o : grid.cpp grid.h util.h data.h eos.h

init.o : init.cpp eos.h grid.h init.h util.h data.h glauber.h random.h

int.o : int.cpp util.h 

reconst.o : reconst.cpp data.h eos.h grid.h reconst.h util.h 

util.o : util.cpp util.h 

glauber.o : glauber.cpp util.h glauber.h data.h random.h

freeze.o : freeze.cpp data.h eos.h grid.h util.h freeze.h int.h

freeze_pseudo.o : freeze_pseudo.cpp data.h eos.h grid.h util.h freeze.h int.h

main.o : main.cpp data.h eos.h evolve.h grid.h init.h util.h glauber.h random.h reconst.h advance.h u_derivative.h

clean:
	rm -f $(OBJS)

#tar:
#	tar cf - makefile eos.cpp evolve.cpp grid.cpp init.cpp main.cpp reconst.cpp util.cpp glauber.cpp eos.h evolve.h grid.h init.h util.h glauber.h data.h freeze.cpp freeze.h freeze_pseudo.cpp input > $(COMMAND).tar

#compress:	tar
#	compress $(COMMAND).tar
