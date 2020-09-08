
MUSIC
--------------------------------------

Compilation
--------------------------------------

The MUSIC code can be compiled using standard `cmake`. 

Alternatively, you can compile the MUSIC code using a makefile. In this way,
please make sure the information about the openMP and the directory of the
GSL library is correct in the `src/GNUmakefile`. Then one can compile the code
by typing `make`.



Run MUSIC with an input file
--------------------------------------

An input file is required that contains the line `EndOfData`, preceded by a
list of parameter names and values, one per line, with parameter names and
values separated by a space or tab. If omitted, each parameter will be assigned
a default value. Example input files can be found under the folder `example_inputfiles/`.

The list of possible parameters and default values, along
with a brief description, can be found in `utilities/music_parameters_dict.py`.

To run hydro simulation, one can start by running the python script
`generate_music_inputfile.py` (can be found in `utilites` folder) to generate input files and job running script.
The job script can run under bash and qsub systems.  For some help information
about `generate_music_inputfile.py`,  one can simply type:

    generate_music_inputfile.py -h



Run MUSIC on multiple CPU cores
--------------------------------------
MUSIC uses openMP for parallelization.  For example, to run on two processors 
and use the sample input file, type:

    export OMP_NUM_THREADS=2
    ./MUSIChydro input_example
Once the prerequisites are installed, you can build the package using:

    make -j 10 #Adjust 10 to the number of cores available.

The result will be an executable named **`MUSIChydro`**.
