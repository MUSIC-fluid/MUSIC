#!/usr/bin/env python

import sys
from os import path

def generate_script(folder_name):
    queue_name = 'sw'
    walltime = '3:00:00'
    ppn = 1
    working_folder = path.join(path.abspath('./'), folder_name)

    script = open(path.join(working_folder, "submit_resonance_job.pbs"), "w")
    if 'y' in folder_name:
        script.write(
"""#!/usr/bin/env bash
#PBS -N %s
#PBS -l walltime=%s
#PBS -l nodes=1:ppn=%d
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -A cqn-654-ad
#PBS -q %s
#PBS -d %s

module add ifort_icc/14.0.4

mpirun -np 1 ./mpihydro music_input_4_y 1>mode_4.log 2>mode_4.err
mpirun -np 1 ./mpihydro music_input_13_y 1>mode_13.log 2>mode_13.err
mpirun -np 1 ./mpihydro music_input_14_y 1>mode_14.log 2>mode_14.err

""" % (folder_name, walltime, ppn, queue_name, working_folder))
    else:
        script.write(
"""#!/usr/bin/env bash
#PBS -N %s
#PBS -l walltime=%s
#PBS -l nodes=1:ppn=%d
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -A cqn-654-ad
#PBS -q %s
#PBS -d %s

module add ifort_icc/14.0.4

mpirun -np 1 ./mpihydro music_input_4 1>mode_4.log 2>mode_4.err
mpirun -np 1 ./mpihydro music_input_13 1>mode_13.log 2>mode_13.err
mpirun -np 1 ./mpihydro music_input_14 1>mode_14.log 2>mode_14.err

""" % (folder_name, walltime, ppn, queue_name, working_folder))

    script.close()


if __name__ == "__main__":
    folder_name = str(sys.argv[1])
    generate_script(folder_name)

