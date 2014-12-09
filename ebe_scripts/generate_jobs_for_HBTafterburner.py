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
      + "\n>>>>> Welcome to the event generator! <<<<<\n" + "-"*80 + normal)

walltime = '1:00:00'
node_folder_name = "node"
code_folder_name = 'HBT_correlation'

file_folder_list = glob(path.join(source_folder, '*'))

for ifolder in range(len(file_folder_list)):
    source_folder_path = file_folder_list[ifolder]
    source_folder_name = source_folder_path.split('/')[-1]
    target_folder = path.join(working_folder, source_folder_name)
    copytree(node_folder_name, target_folder)
    open(path.join(target_folder, "job-%d.pbs" % (ifolder+1)), "w").write(
"""
#!/usr/bin/env bash
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

module add ifort_icc/15.0
(
    cd %s
    mkdir results
    cp %s/tmp/surface.dat results/surface.dat
    cp %s/tmp/input results/input 
    ./HBT.e n_KT=2 KT_min=0.0 KT_max=0.2 >> output.log
    ./HBT.e n_KT=2 KT_min=0.4 KT_max=0.6 >> output.log
    ./HBT.e n_KT=2 KT_min=0.8 KT_max=1.0 >> output.log
    mv results ../HBT_results
)
""" % (ifolder+1, walltime, target_folder, code_folder_name, 
       source_folder_path, source_folder_path)
    )
