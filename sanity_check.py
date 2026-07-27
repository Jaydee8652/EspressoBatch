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
from utils.generic_utils import printToLog as pl, createDirectory as cd, cellVolume, writeCSV
from utils.git_utils import getLocation
from utils.params import *

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)

#Params - can be changed

ecutwfc = 55.0
ecutrho = 440.0
conv_thr = "1.D-6"

batchTarget = 100

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

#Make sure there is a directory to process
cifs_path = os.path.join(homeDirectory, "cifs")
createDirectory(cifs_path, "# WARN - No directory found for .cifs to process.", False)

validated = os.path.join(cifs_path, "validated")
createDirectory(validated, "# WARN - No directory found for .cifs to process. Place .cif files in or replace the newly created directory at ["+validated+"]", True)

cifs = {os.path.splitext(file)[0].replace(".cif", ""): file for file in sorted(os.listdir(validated)) if file.endswith('.cif') and os.path.isfile(os.path.join(validated, file))}

if len(cifs) == 0:#Make sure there are .cifs in the directory
    printToLog("# WARN - No .cif files found to process. Place .cif files in ["+ validated + "]")
    quit()
else:
    printToLog("# INFO - " + str(len(cifs)) + " .cif files found at ["+ validated + "]")
    printToLog("# INFO - Following .cif files are available ["+str(list(cifs.values()))+"]")

#Make sure there is a directory for PSEUDOS
pseudosPath = os.path.join(homeDirectory, "PSEUDOS/")
createDirectory(pseudosPath, "# WARN - No directory found for PSEUDOS. Place .UPF files in or replace the newly created directory at", True)

pseudos = [file for file in os.listdir(pseudosPath) if file.endswith('.UPF') and os.path.isfile(os.path.join(pseudosPath, file))]#Get .UPFs from directory

if len(pseudos) == 0:#Make sure there are .UPFs in the directory
    printToLog("# WARN - No .UPF files found. Place .UPF files in ["+ pseudosPath + "]")
    quit()
else:
    printToLog("# INFO - " + str(len(pseudos)) + " .UPF files found at ["+ pseudosPath + "]")
    psuedo_elements = [file.split(".")[0] for file in pseudos]
    printToLog("# INFO - Following elements accounted for ["+str(list(set(psuedo_elements)))+"]")

#Make sure there is a directory for the generated input files
input_path = os.path.join(homeDirectory, "Sanity_Input_Files")
createDirectory(input_path, "# INFO - Directory for input directories created at ["+ input_path + "]", False)

# Create the SUB submission script
post = os.path.join(os.path.join(homeDirectory,"utils"), "extract_energy.py")

SANITY_SUB = os.path.join(input_path, f"SANITY_SUB")    
with open(SANITY_SUB, "w") as file:
    content = f"""
#!/bin/bash

#SBATCH --job-name=SANITY_CHECK
#SBATCH --mail-type={param_slurmVerbosity}
#SBATCH --mail-user={param_email}
#SBATCH --account={param_account}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=00-23:59
#SBATCH --mem-per-cpu=60G

echo "Running sanity checks [$caselist]"

for case in $caselist
do

    echo "Compound [$case] Starting sanity check"
    cd $case
    srun --cpus-per-task=$SLURM_CPUS_PER_TASK pw.x < $case.in > $case.out
    srun --cpus-per-task=1 --ntasks=1 python3 {post}

    cd ..    
    echo "Compound [$case] Finished sanity check"

done
"""    
    print(content.lstrip("\n"), file=file)
    printToLog("# INFO - Created SANITY_SUB file at ["+SANITY_SUB+"]")

#Create local sheet
localSheet = os.path.join(input_path, "sanity_sheet.csv")
if not os.path.isfile(localSheet):
    printToLog("# INFO - No local .csv found. Attempting to initialise")
    with open(localSheet, 'a') as file:
        file.write("[REFCODE]")

existing_directories = [directory for directory in sorted(os.listdir(input_path)) if os.path.isdir(os.path.join(input_path, directory)) and not directory.startswith(".")]

printToLog("# INFO - Following directories are available to run ["+str(list(existing_directories))+"]")
printToLog("# INFO - Enter integer(s) with spaces between entries ('1 2 3') to choose processes to perform.")
options = {
    "1": "Run cif2cell to produce .in files",
    "2": "Batch ["+str(batchTarget)+"] new sanity check calculations to slurm",
    "0": "All in sequence",
}

for key, value in options.items():
    printToLog(f"# INFO -    [{key}] {value}")
choices = input(">")
invalidInputs = []
regex = re.compile('[^0-9 ]')
choices = regex.sub('', choices).strip().split(" ")
if choices.__contains__("0"):
    choices = list(options)
    choices.remove("0")
    
choices = list(set(choices))
for choice in choices:    
    if not options.__contains__(choice):
        invalidInputs.append(choice)
if len(invalidInputs) > 0:
    printToLog("# WARN - The following inputs ["+str(list(set(invalidInputs)))+"] are not supported")
    quit()
printToLog("# INFO - The following processes have been selected ["+str(sorted(choices,key=int))+"]")

if choices.__contains__("1"):
    for refcode, filename in cifs.items():
        df = pd.read_csv(localSheet)

        printToLog("# INFO - Compound [" + refcode + "] Processing .cif file")
        if existing_directories.__contains__(refcode):
            printToLog("# INFO - Compound ["+refcode+"] Previously processed")
            continue
    
        # Create the directory if it doesn't already exist
        refcodeDirectory = os.path.join(input_path, refcode)
        createDirectory(refcodeDirectory, "# INFO - Compound [" + refcode + "] No directory found, created at", False)
        shutil.copyfile(os.path.join(validated, filename), os.path.join(refcodeDirectory, filename))
    
        # Run cif2cell file generation
        in_path = os.path.join(refcodeDirectory, refcode+".in")
        cif_path = os.path.join(validated, filename)
    
        try:
            subprocess.run(f"cif2cell -f {cif_path} -p quantum-espresso --setup-all -o {in_path}",shell=True)
            printToLog("# INFO - Compound [" + refcode + "] Sucessfully ran cif2cell")
        except subprocess.CalledProcessError as e:
            printToLog("# WARN - Compound [" + refcode + "] Error running cif2cell")
            printToLog(str(e))
            continue
    
        if os.path.exists(in_path):
            with open(in_path) as file:
                lines = file.readlines()
                offset = 0
                for number, line in enumerate(lines.copy(), 0): 
                    number += offset
        
                    lineAtom = line.lstrip().split(" ")[0]
                    if psuedo_elements.__contains__(lineAtom) and len(line.split()) > 3:
                        lines[number] = lines[number].rstrip() + " 0 0 0\n"            
                    if "&SYSTEM" in line:
                        content = f"""
&CONTROL
  calculation = 'scf'
  prefix = '{refcode}'
  outdir = '{refcodeDirectory}/'
  pseudo_dir = '{pseudosPath}'
  !nstep = 0
/

"""
                        lines.insert(number, content.lstrip("\n"))
                        offset += 1
                    if "ntyp" in line:
                        del lines[number+1]
                        content = f"""
  ecutwfc = {ecutwfc}
  ecutrho = {ecutrho}
  vdw_corr = 'grimme-d3'
  dftd3_version = 6
/

&ELECTRONS
  conv_thr = {conv_thr}
/

&IONS
/

"""
                        lines.insert(number+1, content.lstrip("\n"))
                    if "ATOMIC_SPECIES" in line:
                        lines.insert(number,"\n")
                        offset += 1
                    if "ATOMIC_POSITIONS" in line:
                        lines.insert(number,"\n")
                        offset += 1
                    if "K_POINTS" in line:
                        del lines[number+1]
                        lines.insert(number+1,"1 1 1  0 0 0\n")
                    if "_PSEUDO" in line:
                        lineAtom = line.split()[0]
                        if lineAtom+"_PSEUDO" in line:                # Add pseuds
                            if any(lineAtom+'.' in pseudo and 'kjpaw' in pseudo for pseudo in pseudos):
                                for pseudo in pseudos:
                                    if lineAtom+'.' in pseudo and 'kjpaw' in pseudo:
                                        lines[number]=line.replace(lineAtom+"_PSEUDO",pseudo)
                            else:
                                printToLog("# WARN - No .UPF file found for atom ["+str(lineAtom.lstrip())+"], needed by compound ["+refcode+"]")
                                shutil.rmtree(refcodeDirectory)
                                quit()# Deliberate crash if a pseudopotential is not available for all atom types required
                                
            with open(in_path,"w") as file:
                for line in lines:
                    file.write(line)
            printToLog("# INFO - Compound [" + refcode + "] Created .in file at ["+in_path+"]")

            if refcode in df['[REFCODE]'].values:
                printToLog("# INFO - Compound ["+ refcode +"] Already present in sheet")
            else:           
                printToLog("# INFO - Compound ["+ refcode +"] Appending to sheet")
                df = pd.concat([df, pd.DataFrame({"[REFCODE]": [refcode], "[BATCH_started]": " "})], ignore_index=True)
        else:
            printToLog("# WARN - Compound [" + refcode + "] No .in file file found")
            continue
        df.to_csv(localSheet, index=False)

if choices.__contains__("2"):
    df = pd.read_csv(localSheet)  
    
    directories = [directory for directory in sorted(os.listdir(input_path)) if os.path.isdir(os.path.join(input_path, directory)) and directory in df['[REFCODE]'].values]

    df.set_index('[REFCODE]', inplace = True)
    df = df.astype(object)

    unrun = [directory for directory in directories if not str(df.at[directory, "[BATCH_started]"]) == "True" ]    
    unrun = unrun[:min(len(unrun), batchTarget)]
    printToLog(f"# INFO - Attempting to batch [{unrun[0]}->{unrun[-1]}]")
    caselist = ' '.join(unrun)

    if os.path.exists(SANITY_SUB):
        try:
            subprocess.call(f"module load {param_modules}; cd {input_path}; caselist=\"{caselist}\" sbatch SANITY_SUB",shell=True)
            now = str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

            for refcode in unrun: 
                refcodeDirectory = os.path.join(input_path, refcode)
                batch_path = os.path.join(refcodeDirectory, refcode+"_batch.txt")

                with open(batch_path, "a") as batch:
                    writeCSV(df, refcode, "[BATCH_started]", True)
                    writeCSV(df, refcode, "[BATCH_start_time]", now)
                    writeCSV(df, refcode, "[BATCH_location]", getLocation())
                    
                    print("\n# -Batch data\n", file=batch)
                    print("BATCH_started = "+str(True), file=batch)
                    print("BATCH_start_time = "+str(now), file=batch)
                    print("BATCH_location = "+str(getLocation()), file=batch)
            df.to_csv(localSheet)
        except subprocess.CalledProcessError as e:
            printToLog("# WARN - Error batching calculation for compound with refcode ["+refcode+"]")
            printToLog(str(e))
    else:
         printToLog("# WARN - SANITY_SUB not present")
