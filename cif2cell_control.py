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
batch_cap = 16

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

cifs = [os.path.splitext(file)[0].replace(".cif", "") for file in sorted(os.listdir(validated)) if file.endswith('.cif') and os.path.isfile(os.path.join(validated, file))]

if len(cifs) == 0:#Make sure there are .cifs in the directory
    printToLog("# WARN - No .cif files found to process. Place .cif files in ["+ validated + "]")
    quit()
else:
    printToLog("# INFO - " + str(len(cifs)) + " .cif files found at ["+ validated + "]")

#Make sure there is a .csv to read from
qe_params = os.path.join(homeDirectory, "_qe_params.csv")
if not os.path.isfile(qe_params):
    with open(qe_params, 'a') as file:
        print("set_id,test_time,ecutwfc,ecutrho_factor,conv_thr,q_gipaw,calculation,volume_cap,atoms_to_optimise",file=file)
        print("MAIN,10,55.0,8.0,1.D-6,0.01,relax,0,H",file=file)
    printToLog(f"# WARN - No qe_params.csv found. Initialised with default settings at [{qe_params}]")
    quit()

printToLog("# INFO - Enter integer(s) with spaces between entries ('1 2') to choose processes to perform.")
options = {
    "1": "Run test calculations in a slurm job array",
    "2": "Run test calculations in current session",
}

for key, value in options.items():
    printToLog(f"# INFO -    [{key}] {value}")
if len(sys.argv) > 1:
    choices = ' '.join(sys.argv[1:])
else:
    choices = input(">")
 
invalidInputs = []
regex = re.compile('[^0-9 ]')
choices = list(set(regex.sub('', choices).strip().split(" ")))
if len(choices) > 1:
    printToLog("# WARN - Only one process can be run at a time")
    quit()
    
for choice in choices:    
    if not options.__contains__(choice):
        invalidInputs.append(choice)
if len(invalidInputs) > 0:
    printToLog("# WARN - The following inputs ["+str(list(set(invalidInputs)))+"] are not supported")
    quit()
printToLog("# INFO - The following processes have been selected ["+str(sorted(choices,key=int))+"]")


df = pd.read_csv(qe_params, encoding="utf-8-sig")
df.set_index('set_id', inplace = True)

batch_count = 0
for set_id, row in df.iterrows():
    printToLog(f"# INFO - Processing set_id [{set_id}]")

    id_directory = os.path.join(homeDirectory, set_id)
    createDirectory(id_directory, f"# INFO - Directory for paramater set [{set_id}] created at [{id_directory}]", False)
    df.loc[[set_id]].to_csv(os.path.join(id_directory, "local_params.csv"))

    #Make sure there is a directory for the generated input files
    input_files = os.path.join(id_directory, "input_files")
    createDirectory(input_files, f"# INFO - Directory for input directories created at [{input_files}]", False)
    
    unrun = [cif for cif in cifs if not os.path.isdir(os.path.join(input_files, cif))]
    printToLog("# INFO - Following .cif files are available to run ["+str(list(unrun))+"]")

    processing = os.path.join(os.path.join(homeDirectory, "utils"), "qe_cif2cell.py")
    if choices.__contains__("1"):
        seconds = min(grouping,len(unrun)) * float(row['test_time'])
        time_format = time.strftime("00-%H:%M", time.gmtime(seconds))
        
        # Determines sensible grouping to minimise the amount of individual submissions
        if len(unrun) >= grouping:
            grouping = math.ceil((len(unrun) % 100) / math.floor(len(unrun) / grouping)) + grouping
        printToLog(f"# INFO - [{len(unrun)}] .cifs to process for set_id [{set_id}]. Reasonable batch grouping determined to be [{math.ceil(len(unrun) / grouping)}] job(s) each lasting [{time_format}] and containing [{min(grouping,len(unrun))}] compounds")
    
        if round(len(unrun) / grouping) > batch_cap:
              printToLog(f"# WARN - Determined number of job(s) [{round(len(unrun) / grouping)}] is greater than batch cap [{batch_cap}]")     
    
        # Creating slurm submission script
        name = f"_QE_CIF2CELL_{set_id}"
        sub = os.path.join(input_files, name)    
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

echo "Running [{set_id}] qe_cif2cell on [$caselist]"

for case in $caselist
do

    echo "Compound [$case] Starting [{set_id}] qe_cif2cell"
    cd $case
    srun --cpus-per-task=1 --ntasks=1 python3 {processing} prepare_test

    echo "Compound [$case] Running test calculation"
    srun --cpus-per-task=$SLURM_CPUS_PER_TASK pw.x < $case.in > test_$case.out
    srun --cpus-per-task=1 --ntasks=1 python3 {processing}
    
    cd ..    
    echo "Compound [$case] Finished [{set_id}] qe_cif2cell"

done
echo "Finished running [{set_id}] qe_cif2cell on [$caselist]"

"""    
            print(content.lstrip("\n"), file=file)
            printToLog(f"# INFO - Created [{name}] file at [{sub}]")
        
        while len(unrun) > 0 and batch_count <= batch_cap:
            # Get a group of length "grouping" of cifs not yet processed
            unrun = unrun[:min(len(unrun), grouping)]
            printToLog(f"# INFO - [{set_id}] Attempting to batch group of [{min(grouping,len(unrun))}] [{unrun[0]}->{unrun[-1]}]")
            caselist = ' '.join(unrun)
        
            if os.path.exists(sub):
                try:
                    for refcode in unrun: 
                        refcodeDirectory = os.path.join(input_files, refcode)
                        createDirectory(refcodeDirectory, f"# INFO - Compound [{set_id} {refcode}] No directory found, created at", False)
        
                        incomplete = os.path.join(refcodeDirectory, "INCOMPLETE.txt")
                        with open(incomplete, "a") as file:
                            print("WARNING, the presence of this file indicates that the qe_cif2cell process did not run to completion. This input should not be run!", file=file)      
                    subprocess.call(f"module load {param_modules}; cd {input_files}; caselist=\"{caselist}\" sbatch {name}",shell=True)
                    printToLog(f"# INFO - [{set_id}] Successfully batched group of [{min(grouping,len(unrun))}] [{unrun[0]}->{unrun[-1]}]")
        
                    batch_count += 1
                except subprocess.CalledProcessError as e:
                    printToLog("# WARN - Error batching calculation")
                    printToLog(str(e))
            else:
                 printToLog("# WARN - sub not present")
            unrun = [cif for cif in cifs if not os.path.isdir(os.path.join(input_files, cif))]
    elif choices.__contains__("2"):
        for refcode in unrun: 
            printToLog(f"# INFO - Processing compound [{set_id} {refcode}]")
            
            refcodeDirectory = os.path.join(input_files, refcode)
            createDirectory(refcodeDirectory, f"# INFO - Compound [{set_id} {refcode}] No directory found, created at", False)

            incomplete = os.path.join(refcodeDirectory, "INCOMPLETE.txt")
            with open(incomplete, "a") as file:
                print("WARNING, the presence of this file indicates that the qe_cif2cell process did not run to completion. This input should not be run!", file=file)
            try:
                subprocess.call(f"module load {param_modules}; cd {input_files}; cd {refcode}; python3 {processing} prepare_test; pw.x < {refcode}.in > test_{refcode}.out; python3 {processing}",shell=True)
            except subprocess.CalledProcessError as e:
                printToLog(f"# WARN - Error running test calculation for [{refcode}]")
                printToLog(str(e))