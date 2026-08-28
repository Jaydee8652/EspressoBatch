# reprocess.py - Jacob Duddridge

# Runs post_processing.py on all directories containing PWSCF and GIPAW .out files

# All processes are reported to reprocess.log for debugging
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
from utils.generic_utils import printToLog as pl, createDirectory as cd, cellVolume, getQueued
from utils.params import *

ntasks_per_node = 1 #32
mem_per_cpu = 4
grouping = 100
batchCap = 16
runTime = 30 #seconds

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

#Make sure there is a .csv to read from
qe_params = os.path.join(homeDirectory, "qe_params.csv")
if not os.path.isfile(qe_params):
    printToLog("# WARN - No qe_params.csv found.")
    quit()

df = pd.read_csv(qe_params, encoding="utf-8-sig")
df.set_index('set_id', inplace = True)

queued = getQueued(log)

for set_id, row in df.iterrows():
    printToLog(f"# INFO - Processing set_id [{set_id}]")
    
    id_directory = os.path.join(homeDirectory, set_id)
    createDirectory(id_directory, f"# INFO - Directory for paramater set [{set_id}] created at [{id_directory}]", False)
    df.loc[[set_id]].to_csv(os.path.join(id_directory, "local_params.csv"))

    #Make sure there is a directory for the generated input files
    input_files = os.path.join(id_directory, "input_files")
    createDirectory(input_files, f"# INFO - Directory for input directories created at [{input_files}]", False)
    targets = sorted([directory for directory in os.listdir(input_files) if os.path.isdir(os.path.join(input_files, directory)) and not directory.startswith(".") and os.path.isfile(os.path.join(os.path.join(input_files, directory), directory+".out")) and os.path.isfile(os.path.join(os.path.join(input_files, directory), "gipaw."+directory+".out"))])
    
    targets = [refcode for refcode in targets if set_id+"_"+refcode+"_SUB" not in queued]
    
    printToLog("# INFO - Following input directories are available to run ["+str(list(targets))+"]")
    
    time_format = time.strftime("00-%H:%M", time.gmtime(min(grouping,len(targets)) * runTime))
    if len(targets) >= grouping:
        grouping = math.ceil((len(targets) % 100) / math.floor(len(targets) / grouping)) + grouping
    printToLog(f"# INFO - [{len(targets)}] directories to process for set_id [{set_id}]. Reasonable batch grouping determined to be [{math.ceil(len(targets) / grouping)}] job(s) each lasting [{time_format}] and containing [{min(grouping,len(targets))}] compounds")
    
    if round(len(targets) / grouping) > batchCap:
          printToLog(f"# WARN - Determined number of job(s) [{round(len(targets) / grouping)}] is greater than batch cap [{batchCap}]")      
    
    name = f"_REPROCESS_{set_id}"
    sub = os.path.join(input_files, name)    
    processing = os.path.join(os.path.join(homeDirectory, "utils"), "post_processing.py")
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


echo "Rerunning post_processing.py on [$caselist]"

for case in $caselist
do

    echo "Compound [$case] Rerunning post_processing.py"
    cd $case
    srun --cpus-per-task=1 --ntasks=1 python3 {processing}
    cd ..    
    echo "Compound [$case] Finished rerunning post_processing.py"

done

echo "Finished rerunning post_processing.py [$caselist]"
    """    
        print(content.lstrip("\n"), file=file)
        printToLog(f"# INFO - Created {name} file at [{sub}]")
    
    batchCount = 0
    unbatched = targets.copy()
    while len(targets) > 0 and batchCount <= batchCap:
        targets = targets[:min(len(targets), grouping)]
        printToLog(f"# INFO - Attempting to batch group of [{len(targets)}] [{targets[0]}->{targets[-1]}]")
        caselist = ' '.join(targets)
    
        if os.path.exists(sub):
            try:
                subprocess.call(f"module load {param_modules}; cd {input_files}; caselist=\"{caselist}\" sbatch {name}",shell=True)
                printToLog(f"# INFO - Successfully batched group of [{len(targets)}] [{targets[0]}->{targets[-1]}]")
    
                batchCount += 1
            except subprocess.CalledProcessError as e:
                printToLog("# WARN - Error batching calculation")
                printToLog(str(e))
        else:
             printToLog(f"# WARN - {name} sub not present")
        unbatched = [directory for directory in unbatched if directory not in targets]
        targets = unbatched