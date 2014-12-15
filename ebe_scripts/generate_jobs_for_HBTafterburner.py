#! /usr/bin/env python
"""
    This script duplicates the afterburner codes in the node folder and 
    generate a collection of pbs files to be batch-submitted. 
    For efficiency all codes inside node should be compiled.
"""

from sys import argv, exit
from os import makedirs, path, unlink
from shutil import copytree, copy, rmtree
from glob import glob
from numpy import pi

# define colors
purple = "\033[95m"
green = "\033[92m"
blue = "\033[94m"
yellow = "\033[93m"
red = "\033[91m"
normal = "\033[0m"

try:
    source_folder = path.abspath(argv[1])
    working_folder = path.abspath(argv[2])
except(IOError):
    print("Usage: generate_jobs_for_afterburner.py source_folder working_folder")
    exit(1)

print(purple + "\n" + "-"*80 
      + "\n>>>>> Welcome to the event-by-event afterburner generator! <<<<<\n" + "-"*80 + normal)

walltime = '5:00:00'
node_folder_name = "node"
code_folder_name = 'HBT_correlation'
azimuthal_flag = 1
phi_order = 2
module_pi = 2*pi/phi_order

file_folder_list = glob(path.join(source_folder, '*'))

for ifolder in range(len(file_folder_list)):
    source_folder_path = file_folder_list[ifolder]
    source_folder_name = source_folder_path.split('/')[-1]
    target_folder = path.join(working_folder, source_folder_name)
    copytree(node_folder_name, target_folder)
    if azimuthal_flag == 0:
        open(path.join(target_folder, "job-%d.pbs" % (ifolder+1)), "w").write(
"""#!/usr/bin/env bash
#PBS -N after_burner_job-%d
#PBS -l walltime=%s
#PBS -l nodes=1:ppn=2
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -A cqn-654-ac
#PBS -m bea
#PBS -M chunshen1987@gmail.com
#PBS -d %s

module add ifort_icc/14.0.4
(
    cd %s
    mkdir results
    cp %s/tmp/surface.dat results/surface.dat
    cp %s/tmp/input results/input 
    ./HBT.e n_KT=2 KT_min=0.0 KT_max=0.2 azimuthal_flag=0 >> output.log
    ./HBT.e n_KT=2 KT_min=0.4 KT_max=0.6 azimuthal_flag=0 >> output.log
    ./HBT.e n_KT=2 KT_min=0.8 KT_max=1.0 azimuthal_flag=0 >> output.log
    mv results ../HBT_results
)
""" % (ifolder+1, walltime, target_folder, code_folder_name, 
       source_folder_path, source_folder_path)
        )
    else:
        open(path.join(target_folder, "job-%d.pbs" % (ifolder+1)), "w").write(
"""#!/usr/bin/env bash
#PBS -N after_burner_job-%d
#PBS -l walltime=%s
#PBS -l nodes=1:ppn=2
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -A cqn-654-ac
#PBS -m bea
#PBS -M chunshen1987@gmail.com
#PBS -d %s

module add ifort_icc/14.0.4
(
    cd iS
    mkdir results
    cp %s/tmp/surface.dat results/surface.dat
    cp %s/tmp/input results/input 
    ./iS.e
    rand_rotation=$(echo $(($RANDOM %% %d))*%f | bc)
    phi_ev=$(head -n %d results/thermal_211_integrated_vndata.dat | tail -n 1 | awk {'print -atan2($3, $2)/%d'})
    phi_ev=$(echo $phi_ev + $rand_rotation | bc)
    echo $rand_rotation $phi_ev
    cd ../%s
    mv ../iS/results ./
    ./HBT.e n_KT=2 KT_min=0.0 KT_max=0.2 Psi_ev=$phi_ev azimuthal_flag=1 >> output.log
    ./HBT.e n_KT=2 KT_min=0.4 KT_max=0.6 Psi_ev=$phi_ev azimuthal_flag=1 >> output.log
    ./HBT.e n_KT=2 KT_min=0.8 KT_max=1.0 Psi_ev=$phi_ev azimuthal_flag=1 >> output.log
    mv results ../HBT_results
)
""" % (ifolder+1, walltime, target_folder, source_folder_path, 
       source_folder_path, phi_order, module_pi, phi_order+1, phi_order, code_folder_name)
        )
