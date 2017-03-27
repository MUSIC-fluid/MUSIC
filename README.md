#MUSIC               {#mainpage}
===================

MUSIC requires GSL libraries to be installed from
[here](https://www.gnu.org/software/gsl/)
as well as MPI (e.g. MPICH or OpenMPI) and a C++ compiler.  
It only runs properly on POSIX operating systems 
(tested on LINUX and Mac OS X).

## Compile the code

* The MUSIC code can be compiled using standard cmake. 
* Alternatively, you can compile the MUSIC code using a makefile. In this way,
  please make sure the information about the MPI compiler and the directory of
  the GSL library is correct in the src/GNUmakefile. Then one can compile
  the code by typing `make`.

## Run MUSIC with an input file

An input file is required that contains the line "EndOfData", 
preceded by a list of parameter names and values, one per line,
with parameter names and values separated by a space or tab.
If omitted, each parameter will be assigned a default value.  
The list of possible parameters and default values, along
with a brief description, can be found in "music_parameters_dict.py"

To run hydro simulation, one can start by running the python script
"generate_music_inputfile.py" to generate input files and job
running script. The job script can run under bash and qsub systems. 
For some help information about "generate_music_inputfile.py", 
one can simply type:
`generate_music_inputfile.py -h`

### Run MUSIC on multiple CPU cores
MUSIC can be invoked by MPI.  For example, to run on two processors 
and use the sample input file, type:
`mpiexec -n 2 ./mpihydro input_example`

