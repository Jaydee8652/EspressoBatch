# qe_cif2cell.py - Jacob Duddridge

# Generates quantumespresso and gipaw input files from a .cif
# Runs test calculation and updates resource requests accordingly
# Atom fixing is included

# All processes are reported to qe_cif2cell.log for debugging

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

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)

nNodesMax = 192 # Determined by cluster
archive = False # Replace current Input_Files and save existing as an archive file.

ecutwfc = 55.0
ecutrho = 440.0
conv_thr = "1.D-6"
q_gipaw = 0.01

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
inputPath = os.path.join(homeDirectory, "Input_Files")
if os.path.exists(inputPath) and archive:
    now = str(datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
    printToLog("# INFO - Archiving existing file for input directories ["+ inputPath + "]")
    os.rename(inputPath, os.path.join(homeDirectory, "Input_Files_"+now))
createDirectory(inputPath, "# INFO - Directory for input directories created at ["+ inputPath + "]", False)

existing_directories = [directory for directory in os.listdir(inputPath) if os.path.isdir(os.path.join(inputPath, directory)) and not directory.startswith(".")]

printToLog("# INFO - Enter atom types to optimise, with spaces between entries ('H O C'). Enter 'All' to optimise all atoms")
atomsToOptimise = input(">")
invalidInputs = []

# Input sanitisation - user input of atoms to freeze
regex = re.compile('[^a-zA-Z ]')
atomsToOptimise = regex.sub('', atomsToOptimise).strip().split(" ")
atomsToOptimise = list(set([atom.lower().capitalize() for atom in atomsToOptimise]))
if atomsToOptimise.__contains__("All"):
    atomsToOptimise.clear()
    atomsToOptimise.append("All")
for atom in atomsToOptimise:    
    if not psuedo_elements.__contains__(atom) and not atom == "All":
        invalidInputs.append(atom)
if len(invalidInputs) > 0:
    printToLog("# WARN - The following atom types ["+str(list(set(invalidInputs)))+"] are not accounted for by the available .UPF files and have been removed")
atomsToOptimise = [atom for atom in atomsToOptimise if atom not in invalidInputs]
printToLog("# INFO - Following atom types selected to be optimised ["+str(atomsToOptimise)+"]")

for refcode, filename in cifs.items():
    printToLog("# INFO - Compound [" + refcode + "] Processing .cif file")
    if existing_directories.__contains__(refcode):
        printToLog("# INFO - Compound ["+refcode+"] Previously processed")
        continue

    with open(os.path.join(validated, filename), "r") as cif:
        lines = cif.readlines()
        for number, line in enumerate(lines, 0):   
            if not line.startswith("#") or line.startswith(""):
                if "_cell_length_a" in line:
                    cell_a = float(lines[number].split()[1].split("(")[0])
                    cell_b = float(lines[number+1].split()[1].split("(")[0])
                    cell_c = float(lines[number+2].split()[1].split("(")[0])
                    cell_α = math.radians(float(lines[number+3].split()[1].split("(")[0]))
                    cell_β = math.radians(float(lines[number+4].split()[1].split("(")[0]))
                    cell_γ = math.radians(float(lines[number+5].split()[1].split("(")[0]))

                    volume = cellVolume(cell_a, cell_b, cell_c, cell_α, cell_β, cell_γ)
                    printToLog(f"# INFO - Compound [{refcode}] Found to have the following dimensions [a {round(cell_a, 4)}] [b {round(cell_b, 4)}] [c {round(cell_c, 4)}] [α {round(cell_α, 4)}] [β {round(cell_β, 4)}] [γ {round(cell_γ, 4)}] [volume {round(volume, 4)}]")

    # Create the directory if it doesn't already exist
    refcodeDirectory = os.path.join(inputPath, refcode)
    createDirectory(refcodeDirectory, "# INFO - Compound [" + refcode + "] No directory found, created at", False)
    incomplete = os.path.join(refcodeDirectory, "INCOMPLETE.txt")
    with open(incomplete, "a") as file:
        print("WARNING, the presence of this file indicates that the qe_cif2cell process did not run to completion. This input should not be run!", file=file)
    
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
                if not atomsToOptimise.__contains__(lineAtom) and not atomsToOptimise.__contains__("All") and psuedo_elements.__contains__(lineAtom) and len(line.split()) > 3:
                    lines[number] = lines[number].rstrip() + " 0 0 0\n"            
                if "&SYSTEM" in line:
                    content = f"""
&CONTROL
  calculation = 'relax'
  prefix = '{refcode}'
  outdir = '{refcodeDirectory}/'
  pseudo_dir = '{pseudosPath}'
  nstep = 0
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
    
        # Create the .gipaw input file
        gipaw = os.path.join(refcodeDirectory, f"gipaw.{refcode}.in")
        with open(gipaw, "w") as file:
            content = f"""
&inputgipaw
   job = 'nmr'
   prefix = '{refcode}'
   tmp_dir = '{refcodeDirectory}/'
   diagonalization = 'david'
   verbosity = 'high'
   q_gipaw = {q_gipaw}
   spline_ps = .true.
   use_nmr_macroscopic_shape = .true.
/        
"""
            print(content.lstrip(), file=file)
        printToLog("# INFO - Compound [" + refcode + "] Created gipaw .in file at ["+gipaw+"]")
        
        # Run test calculation
        testOutPath = os.path.join(refcodeDirectory, refcode+"_test.out")
        try:
            subprocess.run(f"module load {param_modules}; pw.x < {in_path} > {testOutPath}",shell=True)
            printToLog("# INFO - Compound [" + refcode + "] Successfully ran test command")
        except subprocess.CalledProcessError as e:
            printToLog("# WARN - Compound [" + refcode + "] Error running test command")
            printToLog(str(e))
            continue
            
        # Process test result
        if os.path.exists(testOutPath):
            kPoints = 0
            dynamicalRAM = 0
            
            with open(testOutPath, "r") as test:
                printToLog("# INFO - Compound [" + refcode + "] Processing test output file")
                for line in test:                
                    if len(line) == 0 or line.startswith("#"):
                        continue#   skip blank lines and comments
                    if "number of k points" in line:
                        kPoints = math.ceil(float(re.sub("[^0-9]", "", line)))
                        printToLog("# INFO - Compound [" + refcode + "] Determined to have ["+str(kPoints)+"] k points")    
                    if "Estimated max dynamical RAM per process" in line:
                        dynamicalRAM = (math.ceil(float(re.sub("[^0-9.]", "", line))) * 3)
                        printToLog("# INFO - Compound [" + refcode + "] Determined to use ["+str(dynamicalRAM)+"G] max dynamical RAM per process")    
                        
            if dynamicalRAM > (3.7 * kPoints):
                mult = math.ceil(dynamicalRAM / (3.7 * kPoints))
                kPoints *= mult
                printToLog("# INFO - Compound [" + refcode + "] Has high RAM usage. Multiplied number of tasks per node by ["+str(mult)+"]")
            if kPoints > nNodesMax:
                printToLog("# WARN - Compound [" + refcode + "] Number of tasks per node ["+str(kPoints)+"] above cap of ["+str(nNodesMax)+"]")   
                kPoints = nNodesMax
    
            days = "00"
            if volume > 5000:
                printToLog("# INFO - Compound [" + refcode + "] Volume greater than [5000], 2 extra days allocated.")
                days = "02"
            
            printToLog("# INFO - Compound [" + refcode + "] Number of tasks per node set to ["+str(kPoints)+"]")
            printToLog("# INFO - Compound [" + refcode + "] Memory request set to ["+str(dynamicalRAM)+"G]")
    
            # Create the QE_SUB submission script
            QE_SUB = os.path.join(refcodeDirectory, f"QE_SUB")    
            post = os.path.join(os.path.join(homeDirectory,"utils"), "post_processing.py")
            with open(QE_SUB, "w") as file:
                content = f"""
#!/bin/bash
        
#SBATCH --job-name=[SUB]_{refcode}
#SBATCH --mail-type={param_slurmVerbosity}
#SBATCH --mail-user={param_email}
#SBATCH --account={param_account}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={kPoints}
#SBATCH --cpus-per-task=1
#SBATCH --mem={dynamicalRAM}G
#SBATCH --time={days}-23:59

module load {param_modules}

srun --cpus-per-task=$SLURM_CPUS_PER_TASK pw.x < {refcode}.in > {refcode}.out
srun --cpus-per-task=$SLURM_CPUS_PER_TASK gipaw.x < gipaw.{refcode}.in > gipaw.{refcode}.out
srun --cpus-per-task=1 --ntasks=1 python3 {post}
"""
                print(content.lstrip("\n"),file=file)
            printToLog("# INFO - Compound [" + refcode + "] Created QE_SUB file at ["+QE_SUB+"]")
        
            printToLog("# INFO - Compound [" + refcode + "] Updating .in file")
            with open(in_path) as file:
                lines = file.read().splitlines()
            with open(in_path, "w") as file:
                for line in lines:
                    if "nstep" in line:
                        print(re.sub("  nstep = 0", "  !nstep = 0", line), file=file)
                        printToLog("# INFO - Compound [" + refcode + "] Ready to run real calculation!")
                        os.remove(incomplete)
                    else:
                        print(line, file=file)
        else:
            printToLog("# WARN - Compound [" + refcode + "] No test output found")
            continue
    else:
        printToLog("# WARN - Compound [" + refcode + "] No .in file file found")
        continue
