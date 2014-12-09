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
    working_folder = path.abspath(argv[1])
except(IOError):
    print("Usage: generate_jobs_for_afterburner.py working_folder")
    exit(1)

print(purple + "\n" + "-"*80 
      + "\n>>>>> Welcome to the event generator! <<<<<\n" + "-"*80 + normal)

walltime = '1:00:00'
node_folder_name = "node"
code_folder_name = 'HBT_correlation'

file_folder_list = glob(path.join(working_folder, '*'))

for ifolder in range(len(file_folder_list)):
    target_folder = path.join(file_folder_list[ifolder], 'job-%d' % (ifolder+1))
    copytree(node_folder_name, target_folder)
    open(path.join(target_folder, "job-%d.pbs" % (ifolder+1)), "w").write(
"""
#!/usr/bin/env bash
#PBS -N after_burner_job-%d
#PBS -l walltime=%s
#PBS -l nodes=1:ppn=2
#PBS -j oe
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -A cqn-654-ac
#PBS -d %s
(
    cd %s
    mkdir results
    cp ../../tmp/surface.dat results
    cp ../../tmp/input results
    ./HBT.e n_KT=2 KT_min=0.0 KT_max=0.25
    ./HBT.e n_KT=2 KT_min=0.5 KT_max=0.75
    ./HBT.e n_KT=2 KT_min=1.0 KT_max=1.25
    mv results ../HBT_results
)
""" % (ifolder+1, walltime, target_folder, code_folder_name)
    )