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
from utils.generic_utils import printToLog as pl, createDirectory as cd, getModules

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)

#Params - can be changed
email = "slurmwiddi@gmail.com"

nNodesMax = 192 # Determined by cluster
archive = False # Replace current Input_Files and save existing as an archive file.

#Data For stats
failureCount = 0

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

#Make sure there is a directory to process
validatedPath = os.path.join(os.path.join(homeDirectory, "cifs"), "Validated_CIFs")
createDirectory(validatedPath, "# WARN - No directory found for .cifs to process. Place .cif files in or replace the newly created directory at ["+validatedPath+"]", True)

#Make sure there is a PSEUDOS/ directory for .UPFs
pseudosPath = os.path.join(homeDirectory, "PSEUDOS/")
createDirectory(pseudosPath, "# WARN - No directory found for PSEUDOS. Place .UPF files in or replace the newly created directory at ["+pseudosPath+"]", True)

#Make sure there is a directory for the generated input files
inputPath = os.path.join(homeDirectory, "Input_Files")
if os.path.exists(inputPath) and archive:
    now = str(datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
    printToLog("# INFO - Archiving existing file for input directories ["+ inputPath + "]")
    os.rename(inputPath, os.path.join(homeDirectory, "Input_Files_"+now))
createDirectory(inputPath, "# INFO - Directory for input directories created at ["+ inputPath + "]", False)

existingDirectories = [directory for directory in os.listdir(inputPath) if os.path.isdir(os.path.join(inputPath, directory)) and not directory.startswith(".")]

cifFiles = [file for file in os.listdir(validatedPath) if file.endswith('.cif') and os.path.isfile(os.path.join(validatedPath, file))]#Get .cifs from directory
numberOfCifs = len(cifFiles) # determine number of .cif files
if numberOfCifs == 0:#Make sure there are .cifs in the directory
    printToLog("# WARN - No .cif files found to process. Place .cif files in ["+ validatedPath + "]")
    quit()
else:
    printToLog("# INFO - " + str(numberOfCifs) + " .cif files found at ["+ validatedPath + "]")

pseudosFiles = [file for file in os.listdir(pseudosPath) if file.endswith('.UPF') and os.path.isfile(os.path.join(pseudosPath, file))]#Get .UPFs from directory
numberOfPsueds = len(pseudosFiles) # determine number of .UPF files
psuedoElements = list()
if numberOfPsueds == 0:#Make sure there are .UPFs in the directory
    printToLog("# WARN - No .UPF files found. Place .UPF files in ["+ pseudosPath + "]")
    quit()
else:
    printToLog("# INFO - " + str(numberOfPsueds) + " .UPF files found at ["+ pseudosPath + "]")

    psuedoElements = [file.split(".")[0] for file in pseudosFiles]
    printToLog("# INFO - Following elements accounted for ["+str(list(set(psuedoElements)))+"]")#Atoms with available pseudopotentials

# Get structure data .csv
structureDataPath = os.path.join(homeDirectory, "structure_data.csv")
if os.path.exists(structureDataPath):
    df = pd.read_csv(structureDataPath, encoding="utf-8-sig")
    df.set_index('[REFCODE]', inplace = True)
else:
    printToLog("# WARN - No .csv file found to load compound data. Copy .csv from the CSD to the following path ["+ structureDataPath + "]")
    quit()


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
    if not psuedoElements.__contains__(atom) and not atom == "All":
        invalidInputs.append(atom)
if len(invalidInputs) > 0:
    printToLog("# WARN - The following atom types ["+str(list(set(invalidInputs)))+"] are not accounted for by the available .UPF files and have been removed")
atomsToOptimise = [atom for atom in atomsToOptimise if atom not in invalidInputs]
printToLog("# INFO - Following atom types selected to be optimised ["+str(atomsToOptimise)+"]")

for cif in cifFiles:
    refcode = os.path.splitext(cif)[0]
    printToLog("# INFO - Compound [" + refcode + "] Processing .cif")
    if existingDirectories.__contains__(refcode):
        printToLog("# INFO - Compound ["+refcode+"] Previously processed")
        continue
    
    # Create the directory if it doesn't already exist
    refcodeDirectory = os.path.join(inputPath, refcode)
    createDirectory(refcodeDirectory, "# INFO - Compound [" + refcode + "] No directory found, created at", False)
    incomplete = os.path.join(refcodeDirectory, "INCOMPLETE.txt")
    with open(incomplete, "a") as file:
        print("WARNING, the presence of this file indicates that the qe_cif2cell process did not run to completion. This input should not be run!", file=file)
    
    shutil.copyfile(os.path.join(validatedPath, cif), os.path.join(refcodeDirectory, cif))

    # Create the .gipaw input file
    gipaw = os.path.join(refcodeDirectory, f"gipaw.{refcode}.in")
    with open(gipaw, "w") as f:
        f.write(f"&inputgipaw\n")
        f.write(f"   job = 'nmr'\n")
        f.write(f"   prefix = '{refcode}'\n")
        f.write(f"   tmp_dir = '{refcodeDirectory}/'\n")
        f.write(f"   diagonalization = 'david'\n")
        f.write(f"   verbosity = 'high'\n")
        f.write(f"   q_gipaw = 0.01\n")
        f.write(f"   spline_ps = .true.\n")
        f.write(f"   use_nmr_macroscopic_shape = .true.\n")
        f.write(f"/\n")
    printToLog("# INFO - Compound [" + refcode + "] Created GIPAW input file at ["+gipaw+"]")

    days = "00"
    if float(df.at[refcode, "[_cell_volume]"]) > 5000:
        printToLog("# INFO - Compound [" + refcode + "] Volume greater than [5000], 2 extra days allocated.")
        days = "02"
    
    # Create the QE_SUB submission script
    QE_SUB = os.path.join(refcodeDirectory, f"QE_SUB")        
    with open(QE_SUB, "w") as f:
        f.write(f"#!/bin/bash\n\n")

        #SBATCH --job-name=test
        f.write(f"#SBATCH --job-name={refcode}_SUB\n")
        f.write(f"#SBATCH --mail-type=ALL\n")
        f.write(f"#SBATCH --mail-user={email}\n")
        f.write(f"#SBATCH --account=def-cwiddifi\n")
        f.write(f"#SBATCH --nodes=1\n")
        f.write(f"#SBATCH --ntasks-per-node=24\n")
        f.write(f"#SBATCH --cpus-per-task=1\n")
        f.write(f"#SBATCH --mem=60G\n")
        f.write(f"#SBATCH --time={days}-23:59\n\n")
        f.write(f"module load {getModules()}\n\n")
        f.write(f"srun --cpus-per-task=$SLURM_CPUS_PER_TASK pw.x < {refcode}.in > {refcode}.out\n")
        f.write(f"srun --cpus-per-task=$SLURM_CPUS_PER_TASK gipaw.x < gipaw.{refcode}.in > gipaw.{refcode}.out\n")
    printToLog("# INFO - Compound [" + refcode + "] Created QE_SUB file at ["+QE_SUB+"]")

    # Run cif2cell file generation
    qeInPath = os.path.join(refcodeDirectory, refcode+".in")    
    cif2cell_bash_command = [
        "cif2cell", "-f", os.path.join(validatedPath, cif), "-p", "quantum-espresso", "--setup-all", "-o", qeInPath
    ]
    try:
        subprocess.run(cif2cell_bash_command, check=True)
        printToLog("# INFO - Compound [" + refcode + "] Sucessfully ran cif2cell")
    except subprocess.CalledProcessError as e:
        printToLog("# WARN - Compound [" + refcode + "] Error running cif2cell")
        printToLog(str(e))
        failureCount += 1

    # Freeze atoms not selected to be optimised
    with open(qeInPath) as file:
        lines = file.read().splitlines()
    with open(qeInPath, "w") as file:
        for line in lines:
            lineAtom = line.lstrip().split(" ")[0]

            if not atomsToOptimise.__contains__(lineAtom) and not atomsToOptimise.__contains__("All") and psuedoElements.__contains__(lineAtom) and len(line.split()) > 3:
                print(line + " 0 0 0", file=file)
            else:
                print(line, file=file)

    # .in file setup
    inputFile=open(qeInPath)
    file_content=(inputFile.readlines())
    lineNum=0

    for line in file_content:
        lineNum += 1
        if "&SYSTEM" in line:
            SYSTEM_line = lineNum
        if "ntyp" in line:
            file_content.insert(lineNum,"  dftd3_version = 6,\n")
            file_content.insert(lineNum,"  vdw_corr = 'grimme-d3',\n")
            file_content.insert(lineNum,"  ecutrho = 440.0,\n")
            file_content.insert(lineNum,"  ecutwfc = 55.0,\n")
        if "dftd3_version" in line:
            file_content.insert(lineNum+1,"\n")
            file_content.insert(lineNum+1,"/\n")
            file_content.insert(lineNum+1,"&IONS\n")
            file_content.insert(lineNum+1,"\n")
            file_content.insert(lineNum+1,"/\n")
            file_content.insert(lineNum+1,"  conv_thr = 1.D-6\n")
            file_content.insert(lineNum+1,"&ELECTRONS\n")
            file_content.insert(lineNum+1,"\n")
        if "CELL_PARAMETERS" in line:
            file_content.insert(lineNum+3,"\n")
        if "_PSEUDO" in line:
            lineAtom = line.split()[0]

            # Add pseuds
            if lineAtom+"_PSEUDO" in line:
                if not any(lineAtom+'.' in pseudo and 'kjpaw' in pseudo for pseudo in pseudosFiles):
                    printToLog("# WARN - No .UPF file found for atom ["+str(lineAtom.lstrip())+"], needed by compound ["+refcode+"]")
                    shutil.rmtree(refcodeDirectory)
                    failureCount += 1    
                    quit()# Deliberate crash if a pseudopotential is not available for all atom types required
                else:
                    for pseudo in pseudosFiles:
                        if lineAtom+'.' in pseudo and 'kjpaw' in pseudo:
                            file_content[lineNum-1]=line.replace(lineAtom+"_PSEUDO",pseudo)
            space_line=lineNum

    file_content.insert(space_line,"\n")
    file_content.insert(SYSTEM_line-1,"\n")
    file_content.insert(SYSTEM_line-1,"/\n")
    file_content.insert(SYSTEM_line-1,"  nstep = 0\n")
    file_content.insert(SYSTEM_line-1,"  pseudo_dir = '"+pseudosPath+"',\n")
    file_content.insert(SYSTEM_line-1,"  outdir = '"+refcodeDirectory+"/',\n")
    file_content.insert(SYSTEM_line-1,"  prefix = '"+refcode+"',\n")
    file_content.insert(SYSTEM_line-1,"  calculation = 'relax',\n")
    file_content.insert(SYSTEM_line-1,"&CONTROL\n")

    inputFile.close()

    with open("Modified_input.txt","w") as file:
        for line in file_content:
            file.write(line)
    file.close()

    os.system("mv Modified_input.txt %s" % (qeInPath))
    printToLog("# INFO - Compound [" + refcode + "] Created .in file at ["+qeInPath+"]")

    # Run test calculation
    testOutPath = os.path.join(refcodeDirectory, refcode+"_test.out")
    testCommand = f"module load {getModules()}; pw.x < {qeInPath} > {testOutPath}"
    printToLog("# INFO - Compound [" + refcode + "] Running test command, ["+testCommand+"]")        
    try:
        subprocess.call(testCommand,shell=True)
        printToLog("# INFO - Compound [" + refcode + "] Successfully ran test command")
    except subprocess.CalledProcessError as e:
        printToLog("# WARN - Compound [" + refcode + "] Error running test command")
        printToLog(str(e))
        failureCount += 1

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

        if os.path.exists(QE_SUB):
            printToLog("# INFO - Compound [" + refcode + "] Updating QE_SUB to use values from test calculation")

            post = os.path.join(os.path.join(homeDirectory,"utils"), "post_processing.py")
            with open(QE_SUB) as file:
                lines = file.read().splitlines()
            with open(QE_SUB, "w") as file:
                for line in lines:
                    if "#SBATCH --ntasks-per-node" in line:
                        printToLog("# INFO - Compound [" + refcode + "] Number of tasks per node set to ["+str(kPoints)+"]")
                        print(re.sub("24", str(kPoints), line), file=file)
                    elif "#SBATCH --mem" in line:
                        printToLog("# INFO - Compound [" + refcode + "] Memory request set to ["+str(dynamicalRAM)+"G]")
                        print(re.sub("60", str(dynamicalRAM), line), file=file)
                    else:
                        print(line, file=file)
                file.write(f"srun --cpus-per-task=1 --ntasks=1 python3 {post}")
        else:
            printToLog("# WARN - Compound [" + refcode + "] No QE_SUB file found")
            failureCount += 1

        if os.path.exists(qeInPath):
            printToLog("# INFO - Compound [" + refcode + "] Updating .in file")

            with open(qeInPath) as file:
                lines = file.read().splitlines()
            with open(qeInPath, "w") as file:
                for line in lines:
                    if "nstep" in line:
                        print(re.sub("  nstep = 0", "  !nstep = 0", line), file=file)
                        printToLog("# INFO - Compound [" + refcode + "] Ready to run real calculation!")
                        os.remove(incomplete)
                    else:
                        print(line, file=file)
        else:
            printToLog("# WARN - Compound [" + refcode + "] No .in file file found")
            failureCount += 1
        
    else:
        printToLog("# WARN - Compound [" + refcode + "] No test output found")
        failureCount += 1

printToLog("# INFO - Process complete, ["+str(numberOfCifs-failureCount)+"] .cifs successessfully processed vs ["+str(failureCount)+"] failures, a % success of ["+str(round(((numberOfCifs-failureCount)/numberOfCifs)*100, 3))+"%]")