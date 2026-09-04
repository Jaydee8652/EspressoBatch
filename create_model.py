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

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)

ntasks_per_node = 1 #32
mem_per_cpu = 4
grouping = 100
batchCap = 16

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

#Make sure there is a directory to sort
models = os.path.join(homeDirectory, "models")
createDirectory(models, "# INFO - No directory found for models, created at", False)

printToLog("# INFO - Enter name of model")
if len(sys.argv) > 1:
    model_name = ' '.join(sys.argv[1:])
else:
    model_name = input(">")

if len(model_name.split()) > 1:
    printToLog(f"# WARN - Model name [{model_name}] cannot include spaces")
    quit()

printToLog(f"# INFO - Creating model with name [{model_name}]")

directory = os.path.join(models, model_name)
createDirectory(directory, "# INFO - No directory found for model, created at", False)

name = f"_TRAIN_{model_name}"
sub = os.path.join(directory, name)    
processing = os.path.join(os.path.join(homeDirectory, "utils"), "train_model.py")
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
#SBATCH --time=00-23:59
#SBATCH --mem-per-cpu={mem_per_cpu}G

srun --cpus-per-task=1 --ntasks=1 python3 {processing}
"""    
    print(content.lstrip("\n"), file=file)
    printToLog(f"# INFO - Created {name} file at [{sub}]")
printToLog(f"# INFO - Attempting to batch model training for [{model_name}]")
if os.path.exists(sub):
    try:
        subprocess.call(f"module load {param_modules}; cd {directory}; sbatch {name}",shell=True)
        printToLog(f"# INFO - Successfully batched model training for [{model_name}]")
    except subprocess.CalledProcessError as e:
        printToLog("# WARN - Error batching model training")
        printToLog(str(e))
else:
     printToLog(f"# WARN - {name} sub not present")