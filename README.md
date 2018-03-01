=======
MUSIC
======================================

Compilation
--------------------------------------

The MUSIC code can be compiled using standard `cmake`. 

Alternatively, you can compile the MUSIC code using a makefile. In this way,
please make sure the information about the MPI compiler and the directory of the
GSL library is correct in the `src/GNUmakefile`. Then one can compile the code
by typing `make`.



Run MUSIC with an input file
--------------------------------------

An input file is required that contains the line `EndOfData`, preceded by a
list of parameter names and values, one per line, with parameter names and
values separated by a space or tab. If omitted, each parameter will be assigned
a default value. The list of possible parameters and default values, along
with a brief description, can be found in `music_parameters_dict.py`.

To run hydro simulation, one can start by running the python script
`generate_music_inputfile.py` to generate input files and job running script.
The job script can run under bash and qsub systems.  For some help information
about `generate_music_inputfile.py`,  one can simply type:

    generate_music_inputfile.py -h



Run MUSIC on multiple CPU cores
--------------------------------------
MUSIC can be invoked by MPI.  For example, to run on two processors 
and use the sample input file, type:

    mpiexec -n 2 ./mpihydro input_example
=======
Once the prerequisites are installed, you can build the package using:

    make -j 10 #Adjust 10 to the number of cores available.

The result will be an executable named **`mpihydro.exe`**.



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

