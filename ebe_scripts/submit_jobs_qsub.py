#! /usr/bin/env python
"""
    submit the pbs jobs using qsub
"""

from sys import argv, exit
from os import path, listdir
from glob import glob
from subprocess import call

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

file_folder_list = glob(path.join(working_folder, '*'))

for ifolder in range(len(file_folder_list)):
    target_folder = path.abspath(file_folder_list[ifolder])
    for aFile in listdir(target_folder):
        if path.splitext(aFile)[1].lower() == '.pbs':
            commandString = "qsub %s" % aFile
            print("Submitting %s in %s ..." % (aFile, target_folder))
            call(commandString, shell=True, cwd=target_folder)

print("Job submission done.")
