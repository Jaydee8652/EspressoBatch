# cif2cell_control.py - Jacob Duddridge

# Generates quantumespresso and gipaw input files from a .cif
# Runs test calculation and updates resource requests accordingly
# Atom fixing is included

#Imports
import __main__ as main
import os
import subprocess
import sys
import shutil
import re 
import math
import datetime
import time
import pandas as pd
from utils.generic_utils import printToLog as pl, createDirectory as cd, cellVolume
from utils.params import *

ntasks_per_node = 1
mem_per_cpu = 4
grouping = 1000
batchCap = 16

testTime = 10 #seconds

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

#Make sure there is a directory for PSEUDOS
pseudosPath = os.path.join(os.path.join(os.path.join(homeDirectory, "utils"), "data"), "PSEUDOS/")
createDirectory(pseudosPath, "# WARN - No directory found for PSEUDOS. Place .UPF files in or replace the newly created directory at", True)

pseudos = [file for file in os.listdir(pseudosPath) if file.endswith('.UPF') and os.path.isfile(os.path.join(pseudosPath, file))]#Get .UPFs from directory

if len(pseudos) == 0:#Make sure there are .UPFs in the directory
    printToLog("# WARN - No .UPF files found. Place .UPF files in ["+ pseudosPath + "]")
    quit()
else:
    printToLog("# INFO - " + str(len(pseudos)) + " .UPF files found at ["+ pseudosPath + "]")
    psuedo_elements = [file.split(".")[0] for file in pseudos]
    printToLog("# INFO - Following elements accounted for ["+str(list(set(psuedo_elements)))+"]")

#Make sure there is a directory to process
cifs_path = os.path.join(homeDirectory, "cifs")
createDirectory(cifs_path, "# WARN - No directory found for .cifs to process.", False)

validated = os.path.join(cifs_path, "validated")
createDirectory(validated, "# WARN - No directory found for .cifs to process. Place .cif files in or replace the newly created directory at ["+validated+"]", True)

#Make sure there is a directory for the generated input files
input_path = os.path.join(homeDirectory, "Input_Files")
createDirectory(input_path, "# INFO - Directory for input directories created at ["+ input_path + "]", False)

cifs = [os.path.splitext(file)[0].replace(".cif", "") for file in sorted(os.listdir(validated)) if file.endswith('.cif') and os.path.isfile(os.path.join(validated, file))]

if len(cifs) == 0:#Make sure there are .cifs in the directory
    printToLog("# WARN - No .cif files found to process. Place .cif files in ["+ validated + "]")
    quit()
else:
    printToLog("# INFO - " + str(len(cifs)) + " .cif files found at ["+ validated + "]")

unrun = [cif for cif in cifs if not os.path.isdir(os.path.join(input_path, cif))]
printToLog("# INFO - Following .cif files are available to run ["+str(list(unrun))+"]")

seconds = min(grouping,len(unrun)) * testTime
time_format = time.strftime("00-%H:%M", time.gmtime(seconds))

# Determines sensible grouping to minimise the amount of individual submissions
if len(unrun) >= grouping:
    grouping = math.ceil((len(unrun) % 100) / math.floor(len(unrun) / grouping)) + grouping
printToLog(f"# INFO - [{len(unrun)}] .cifs to process. Reasonable batch grouping determined to be [{math.ceil(len(unrun) / grouping)}] job(s) each lasting [{time_format}] and containing [{min(grouping,len(unrun))}] compounds")

if round(len(unrun) / grouping) > batchCap:
      printToLog(f"# WARN - Determined number of job(s) [{round(len(unrun) / grouping)}] is greater than batch cap [{batchCap}]")      

printToLog("# INFO - Enter atom types to optimise, with spaces between entries ('H O C'). Enter 'All' to optimise all atoms")
if len(sys.argv) > 1:
    atomsToOptimise = ' '.join(sys.argv[1:]) # Can take command line inputs ie "python3 cif2cell_control.py H C O"
else:
    atomsToOptimise = input(">")
invalidInputs = []

# Input sanitisation - user input of atoms to freeze
regex = re.compile('[^a-zA-Z ]')
atomsToOptimise = regex.sub('', atomsToOptimise).strip().split(" ")
atomsToOptimise = list(set([atom.lower().capitalize() for atom in atomsToOptimise]))
if atomsToOptimise == [""]: # If the input is blank add "None", needs to have a value for later
    atomsToOptimise.append("None")

if atomsToOptimise.__contains__("All"): # If there is "All", remove everything else. "H C All" -> "All"
    atomsToOptimise.clear()
    atomsToOptimise.append("All")
if atomsToOptimise.__contains__("None"): # If there is "None", remove everything else. "H C None" -> "None"
    atomsToOptimise.clear()
    atomsToOptimise.append("None")
    
for atom in atomsToOptimise:    
    if not psuedo_elements.__contains__(atom) and not atom == "All" and not atom == "None":
        invalidInputs.append(atom)

# If we don't have the pseudo for the atom, scrap it, we couldn't run it anyway
if len(invalidInputs) > 0:
    printToLog("# WARN - The following atom types ["+str(list(set(invalidInputs)))+"] are not accounted for by the available .UPF files and have been removed")
atomsToOptimise = [atom for atom in atomsToOptimise if atom not in invalidInputs]
printToLog("# INFO - Following atom types selected to be optimised ["+str(atomsToOptimise)+"]")

# Creating slurm submission script
name = f"_QE_CIF2CELL"
sub = os.path.join(input_path, name)    
processing = os.path.join(os.path.join(homeDirectory, "utils"), "qe_cif2cell.py")
with open(sub, "w") as file:
    content = f"""
#!/bin/bash

#SBATCH --job-name={name}
#SBATCH --mail-type={param_slurmVerbosity}
#SBATCH --mail-user={param_email}
#SBATCH --account={param_account}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={ntasks_per_node}
#SBATCH --cpus-per-task=1
#SBATCH --time={time_format}
#SBATCH --mem-per-cpu={mem_per_cpu}G


echo "Running qe_cif2cell checks [$caselist]"

for case in $caselist
do

    echo "Compound [$case] Starting qe_cif2cell"
    cd $case
    srun --cpus-per-task=1 --ntasks=1 python3 {processing} {" ".join(atomsToOptimise)}
    
    echo "Compound [$case] Running test calculation"
    srun --cpus-per-task=$SLURM_CPUS_PER_TASK pw.x < $case.in > test_$case.out
    srun --cpus-per-task=1 --ntasks=1 python3 {processing}
    
    cd ..    
    echo "Compound [$case] Finished qe_cif2cell"

done
"""    
    print(content.lstrip("\n"), file=file)
    printToLog(f"# INFO - Created [{name}] file at [{sub}]")


batchCount = 0
while len(unrun) > 0 and batchCount <= batchCap:
    # Get a group of length "grouping" of cifs not yet processed
    unrun = unrun[:min(len(unrun), grouping)]
    printToLog(f"# INFO - Attempting to batch group of [{min(grouping,len(unrun))}] [{unrun[0]}->{unrun[-1]}]")
    caselist = ' '.join(unrun)

    if os.path.exists(sub):
        try:
            for refcode in unrun: 
                refcodeDirectory = os.path.join(input_path, refcode)
                createDirectory(refcodeDirectory, "# INFO - Compound [" + refcode + "] No directory found, created at", False)

                incomplete = os.path.join(refcodeDirectory, "INCOMPLETE.txt")
                with open(incomplete, "a") as file:
                    print("WARNING, the presence of this file indicates that the qe_cif2cell process did not run to completion. This input should not be run!", file=file)      
            subprocess.call(f"module load {param_modules}; cd {input_path}; caselist=\"{caselist}\" sbatch _QE_CIF2CELL",shell=True)
            #subprocess.call(f"module load {param_modules}; cd {input_path}; caselist=\"{caselist}\" {name}",shell=True)

            printToLog(f"# INFO - Successfully batched group of [{min(grouping,len(unrun))}] [{unrun[0]}->{unrun[-1]}]")

            batchCount += 1
        except subprocess.CalledProcessError as e:
            printToLog("# WARN - Error batching calculation for compound with refcode ["+refcode+"]")
            printToLog(str(e))
    else:
         printToLog("# WARN - sub not present")
    unrun = [cif for cif in cifs if not os.path.isdir(os.path.join(input_path, cif))]