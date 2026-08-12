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
ntasks_per_node = 16 #32
mem_per_cpu = 60
ecutwfc = 30.0 #55.0
ecutrho = 240.0 #440.0 
conv_thr = "1.D-6"

grouping = 1000
batchCap = 16

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

#Make sure there is a directory to process
cifs_path = os.path.join(homeDirectory, "cifs")
createDirectory(cifs_path, "# WARN - No directory found for .cifs to process.", False)

validated = os.path.join(cifs_path, "validated")
createDirectory(validated, "# WARN - No directory found for .cifs to process. Place .cif files in or replace the newly created directory at", True)
validated_count = len(os.listdir(validated))


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
createDirectory(input_path, "# INFO - Directory for input directories created at ", False)

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
#SBATCH --ntasks-per-node={ntasks_per_node}
#SBATCH --cpus-per-task=1
#SBATCH --time=00-23:59
#SBATCH --mem-per-cpu={mem_per_cpu}G

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

df = pd.read_csv(localSheet)  
directories = [directory for directory in sorted(os.listdir(input_path)) if os.path.isdir(os.path.join(input_path, directory)) and directory in df['[REFCODE]'].values]

df.set_index('[REFCODE]', inplace = True)
df = df.astype(object)
unrun = [directory for directory in directories if not str(df.at[directory, "[BATCH_started]"]) == "True" ]    

printToLog("# INFO - Following directories are available to run ["+str(list(unrun))+"]")
if len(unrun) >= grouping:
    grouping = math.ceil((len(unrun) % 100) / round(len(unrun) / grouping)) + grouping
printToLog(f"# INFO - [{len(unrun)}] checks to run. Reasonable batch grouping determined to be [{round(len(unrun) / grouping)}] job(s) each containing [{grouping}] calculations")

if round(len(unrun) / grouping) > batchCap:
      printToLog(f"# WARN - Determined number of job(s) [{round(len(unrun) / grouping)}] is greater than batch cap [{batchCap}]")      

printToLog("# INFO - Enter integer to choose process to perform.")
options = {
    "1": f"CALC NODE | Run cif2cell to produce .in files",
    "2": f"HEAD NODE | Batch [{len(unrun)}] ([{round(len(unrun) / grouping)}] jobs each containing [{grouping}]) sanity check calculations to slurm",
    "3": f"Calculate relative energies for all outputs in kJ mol⁻¹ molecule⁻¹",
    "4": f"Discard .cifs from 'cifs/validated' not marked in 'Sanity_Input_Files/sanity_sheet.csv' with [validated] = 'True'",
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

choices = choices[0]
for choice in choices: 
    if not options.__contains__(choice):
        invalidInputs.append(choice)
if len(invalidInputs) > 0:
    printToLog("# WARN - The following input ["+str(list(set(invalidInputs)))+"] is not supported")
    quit()
printToLog("# INFO - The following process has been selected ["+str(sorted(choices,key=int))+"]")

if choices.__contains__("1"):
    for refcode, filename in cifs.items():
        df = pd.read_csv(localSheet)

        printToLog("# INFO - Compound [" + refcode + "] Processing .cif file")
        if existing_directories.__contains__(refcode):
            printToLog("# INFO - Compound ["+refcode+"] Previously processed")
            if not (refcode in df['[REFCODE]'].values):
                df = pd.concat([df, pd.DataFrame({"[REFCODE]": [refcode], "[BATCH_started]": " "})], ignore_index=True)
                df.to_csv(localSheet, index=False)
                print("# WARN - Compound ["+refcode+"] Not in sheet. Appending")
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
    batchCount = 0
    while len(unrun) > 0 and batchCount <= batchCap:
        unrun = unrun[:min(len(unrun), grouping)]
        printToLog(f"# INFO - Attempting to batch group of [{grouping}] [{unrun[0]}->{unrun[-1]}]")
        caselist = ' '.join(unrun)
    
        if os.path.exists(SANITY_SUB):
            try:
                subprocess.call(f"module load {param_modules}; cd {input_path}; caselist=\"{caselist}\" sbatch SANITY_SUB",shell=True)
                now = str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                batchCount += 1

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
        unrun = [directory for directory in directories if not str(df.at[directory, "[BATCH_started]"]) == "True" ]

if choices.__contains__("3"):
    minimmum_energy = float(df['[SCF_final_energy_kJ_mol⁻¹_molecule⁻¹]'].min())
    printToLog(f"# INFO - Minimum energy determined to be [{minimmum_energy} kJmol⁻¹ molecule⁻¹] [{df['[SCF_final_energy_kJ_mol⁻¹_molecule⁻¹]'].idxmin()}]")

    for refcode, row in df.iterrows():
        relative_energy = minimmum_energy - float(row['[SCF_final_energy_kJ_mol⁻¹_molecule⁻¹]'])
        printToLog(f"# INFO - Compound [{refcode}] Relative energy determined to be [{relative_energy} kJ mol⁻¹ molecule⁻¹]")
        writeCSV(df, refcode, "[relative_energy_kJ_mol⁻¹_molecule⁻¹]", relative_energy)
        df.to_csv(localSheet)
    df["[validated]"] = ""
    df.to_csv(localSheet)

if choices.__contains__("4"):
    high_energy = os.path.join(cifs_path, "high_energy")
    createDirectory(high_energy, "# WARN - No directory found for high energy .cifs, created at", False)
    high_energy_count = len(os.listdir(high_energy))

    time = str(datetime.datetime.now().strftime("[%Y-%m-%d_%H-%M-%S]"))
    validated_backup = os.path.join(cifs_path, f"validated_backup_{time}")
    shutil.copytree(validated, validated_backup, dirs_exist_ok=True)

    for refcode, row in df.iterrows():
        if not str(row['[validated]']).lower().capitalize() == "True":
            shutil.move(os.path.join(validated, f"{refcode}.cif"), os.path.join(high_energy, f"{refcode}.cif"))
            printToLog(f"# INFO - Compound [{refcode}] Discarded, not flagged to be retained")

    printToLog(f"# INFO - [{'{: >3}'.format(len(os.listdir(high_energy)) - high_energy_count)}] compounds discarded")
    
    surviving = len(os.listdir(validated))
    printToLog(f"# INFO - Sort complete, [{surviving}] of [{validated_count}] compounds validated, with a % survival of [{round((surviving / len(cifs)) * 100, 3)} %]")
