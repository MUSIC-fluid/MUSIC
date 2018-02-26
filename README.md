MUSIC
======================================

Compilation
--------------------------------------

MUSIC requires [GSL libraries](https://www.gnu.org/software/gsl/) to be
installed as well as MPI (e.g. MPICH or OpenMPI) and a C++ compiler. It only
runs properly on POSIX operating systems (tested on LINUX and Mac OS X).

The following command installs the necessary packages on Ubuntu 17.10, and
possible other systems:

    sudo apt install libopenmpi-dev libgsl-dev

Once the prerequisites are installed, you can build the package using:

    make -j 10 #Adjust 10 to the number of cores available.

The result will be an executable named `mpihydro`.



Input
--------------------------------------


An input file is required that contains the line "EndOfData", 
preceded by a list of parameter names and values, one per line,
with parameter names and values separated by a space or tab.
If omitted, each parameter will be assigned a default value.  
The list of possible parameters and default values, along
with a brief description, can be found by reading the function 
"ReadInData2" in main.cpp.
A simple example is included in file "input_sample", which runs
a quick hydro event on a coarse grid, and can be used to verify 
that the code compiles and runs without errors.

MUSIC can be invoked by MPI.  For example, to run on two processors 
and use the sample input file, type:

    mpiexec -n 2 ./mpihydro input_sample

