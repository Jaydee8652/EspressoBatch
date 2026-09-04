# post_processing.py - Jacob Duddridge

# Extracts desired data from PWSCF and gipaw .out files to a summary file, saved here and in a dedicated directory
# Automatically removes symmetry equivalent atoms from the gipaw output

# Establishes connectivity relationships and produces .mol2 files for a single cell and a single molecule

# Automatically run by slurm within the input directory

import os
import sys
import shutil
import re 
import math
import pandas as pd
import io
import csv
import datetime
import time          
import subprocess
import numpy as np

#Jank thing to fix the path. Very annoying artefact of running python scripts by absolute path.
sys.path[0] = sys.path[0][:-6] + sys.path[0][-6:].replace("/utils", "")

from utils.generic_utils import printToLog as pl, createDirectory as cd, removeDirectory as rd, writeCSV, isQueued, parseAlat, cellVolume, mol2Creator
from utils.params import *
#from generic_utils import printToLog as pl, createDirectory as cd, removeDirectory as rd, writeCSV, isQueued, parseAlat, cellVolume, mol2Creator
#from utils.params import *

#Params - can be modified
tolerance = 0.01

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)    
def createDirectory(path, text, exit):
    cd(log, path, text, exit)
def removeDirectory(path, text):
    rd(log, path, text)

CIF_symmetryElements = []

PWSCF_ecutwfc = ""
PWSCF_ecutrho = "" 
PWSCF_conv_thr = ""
PWSCF_version = ""
PWSCF_start_time = ""
PWSCF_end_time = ""
PWSCF_numberMPI = ""
PWSCF_numberThreads = ""
PWSCF_RG = ""
PWSCF_estimatedRAM = ""
PWSCF_scfCycles = ""
PWSCF_bfgsSteps = ""
PWSCF_finalEnergy = ""
PWSCF_CPU_time = ""
PWSCF_WALL_time = ""

PWSCF_done = False

GIPAW_q_gipaw = ""
GIPAW_version = ""
GIPAW_start_time = ""
GIPAW_end_time = ""
GIPAW_numberMPI = ""
GIPAW_numberThreads = ""
GIPAW_RG = ""
GIPAW_mscPPM = ""
GIPAW_msCorrection = []
GIPAW_CPU_time = ""
GIPAW_WALL_time = ""

GIPAW_done = False

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
refcodeDirectory = os.getcwd()#Directory where we are
set_directory = os.path.split(os.path.split(refcodeDirectory)[0])[0]
homeDirectory = os.path.split(set_directory)[0]

logs = os.path.join(refcodeDirectory, "logs")
if not os.path.exists(logs):
    os.makedirs(logs)

refcode = os.path.basename(refcodeDirectory)
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Compound ["+refcode+"] Starting "+str(os.path.basename(sys.argv[0]).split(".")[0])+" in ["+str(refcodeDirectory)+"]")

local_params = os.path.join(set_directory, "local_params.csv")
df = pd.read_csv(local_params, encoding="utf-8-sig")

#params
set_id = str(df['set_id'].iloc[0])

utils = os.path.join(homeDirectory,"utils")
atom_data = os.path.join(os.path.join(utils, "data"), "atom_data.csv")
if os.path.exists(atom_data):
    df = pd.read_csv(atom_data, encoding="utf-8-sig")
    df.set_index('Symbol', inplace = True)
else:
    printToLog("# WARN - No .csv file found to load atom data.")
    quit()

output_files = os.path.join(set_directory, "output_files")
createDirectory(output_files, "# INFO - No directory found for output files, created at", False)

summary_files = os.path.join(set_directory, "summary_files")
createDirectory(summary_files, "# INFO - No directory found for summary files, created at ", False)

out = os.path.join(refcodeDirectory, refcode+".out")    
cif = os.path.join(refcodeDirectory, refcode+".cif")
opt_cif = os.path.join(refcodeDirectory, refcode+"_opt.cif")
super_mol2 = os.path.join(refcodeDirectory, refcode+"_super.mol2")
summaryPath = os.path.join(refcodeDirectory, refcode+"_summary.txt")
summaryCopyPath =os.path.join(summary_files, refcode+"_summary.txt")

cell_mol2 = os.path.join(refcodeDirectory, refcode+"_cell.mol2")
molecule_mol2 = os.path.join(refcodeDirectory, refcode+"_molecule.mol2")

if os.path.exists(opt_cif):
    printToLog("# INFO - Compound ["+ refcode +"] Cleaning existing optimised _super.cif file ["+ opt_cif + "]")
    os.remove(opt_cif)
if os.path.exists(super_mol2):
    printToLog("# INFO - Compound ["+ refcode +"] Cleaning existing _super.mol2 file ["+ super_mol2 + "]")
    os.remove(super_mol2)
if os.path.exists(cell_mol2):
    printToLog("# INFO - Compound ["+ refcode +"] Cleaning existing _cell.mol2 file ["+ cell_mol2 + "]")
    os.remove(cell_mol2)
if os.path.exists(molecule_mol2):
    printToLog("# INFO - Compound ["+ refcode +"] Cleaning existing _molecule.mol2 file ["+ molecule_mol2 + "]")
    os.remove(molecule_mol2)
if os.path.exists(summaryPath):
    printToLog("# INFO - Compound ["+ refcode +"] Cleaning existing summary file ["+ summaryPath + "]")
    os.remove(summaryPath)
if os.path.exists(summaryCopyPath):
    printToLog("# INFO - Compound ["+ refcode +"] Cleaning copied summary file ["+ summaryCopyPath + "]")
    os.remove(summaryCopyPath)

# Remove .save directory
saveDirectory = os.path.join(refcodeDirectory, refcode+".save")
if os.path.exists(saveDirectory):
    printToLog("# INFO - Compound ["+ refcode +"] Removing .save file ["+ saveDirectory + "]")
    shutil.rmtree(saveDirectory)

# Remove .wfc files
regex = re.compile('[^a-zA-Z.]')
wfcFiles = [file for file in os.listdir(refcodeDirectory) if regex.sub('', file).endswith('.wfc') and os.path.isfile(os.path.join(refcodeDirectory, file))]
for wfc in wfcFiles:
    os.remove(os.path.join(refcodeDirectory, wfc))
printToLog(f"# INFO - Compound [{refcode}] Removing [{len(wfcFiles)}] .wfc files")

printToLog("# INFO - Compound ["+ refcode +"] Populating summary file ["+ summaryPath + "]")
with open(summaryPath, "a") as summary:
    print("#Output summary for compound with refcode ["+refcode+"]", file=summary)
    
    #REFCODE_batch.txt
    batch = os.path.join(refcodeDirectory, refcode+"_batch.txt")
    if os.path.isfile(batch):
        now = str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

        with open(batch) as file:
            read = file.read()
        with open(batch, "a") as file:            
            if not read.__contains__("BATCH_done") and isQueued(log, set_id+"_"+refcode):
                print("BATCH_end_time = "+str(now), file=file)       
                print("BATCH_done = "+str(True), file=file)            
        with open(batch) as file:
            lines = file.read().splitlines()
            for line in lines:
                print(line, file=summary)
    else:
        print("# WARN - No batch file found for compound with refcode ["+refcode+"]", file=summary)
        printToLog("# WARN - Compound ["+refcode+"] No batch file found")
        quit()

    # Get QE_SUB
    QE_SUB = os.path.join(refcodeDirectory,"QE_SUB")
    if os.path.isfile(QE_SUB):
        with open(QE_SUB) as file:
            lines = file.read().splitlines()
            for line in lines:
                value = line[line.find("=")+1:].strip()

                if "#SBATCH --nodes" in line:
                    print(f"BATCH_nodes = {value}", file=summary)

                elif "#SBATCH --ntasks-per-node" in line:
                    print(f"BATCH_ntasks_per_node = {value}", file=summary)

                elif "#SBATCH --cpus-per-task" in line:
                    print(f"BATCH_cpus_per_task = {value}", file=summary)

                elif "#SBATCH --mem" in line:
                    print(f"BATCH_memory = {value}", file=summary)

                elif "#SBATCH --time" in line:
                    print(f"BATCH_alloted_time = {value}", file=summary)
    else:
        print("# WARN - No QE_SUB file found for compound with refcode ["+refcode+"]", file=summary)
        printToLog("# WARN - Compound ["+refcode+"] No QE_SUB file found")
        quit()
        
    # Get .cif
    if os.path.isfile(cif):
        with open(cif) as file:
            lines = file.readlines()
            if "_symmetry_equiv_pos_site_id\n" in lines:
                for line in lines: 
                    if line[0].isdigit():
                        CIF_symmetryElements.append(line.split()[1])  
            else:
                sub = ''.join(lines).split('_symmetry_equiv_pos_as_xyz\n')[1].split('\nloop_')[0].strip().split("\n")
                for number, line in enumerate(sub):
                    CIF_symmetryElements.append((line.strip("'")))
    
            printToLog(f"# INFO - Compound [{refcode}] Has the following symmetry operations [{CIF_symmetryElements}]")  
    else:
        print("# WARN - No .cif file found for compound with refcode ["+refcode+"]", file=summary)
        printToLog("# WARN - Compound ["+refcode+"] No .cif file found")
        quit()
        
    # Get .in
    pwscfIn = os.path.join(refcodeDirectory, refcode+".in")
    if os.path.isfile(pwscfIn):
        print("\n# -PWSCF params-\n", file=summary)            
        with open(pwscfIn, "r") as file:
            reduction_factor = 1
            
            start = 0
            end = 0
    
            lines = file.readlines()
            for number, line in enumerate(lines, 1):  
                if "&SYSTEM" in line:
                    Alat = float(re.sub("[^0-9.]", "", lines[number+1]).strip())
                elif "ecutwfc" in line:
                    PWSCF_ecutwfc = float(re.sub("[^0-9.]", "", line).strip())

                    print("PWSCF_ecutwfc = "+str(PWSCF_ecutwfc), file=summary)
                elif "ecutrho" in line:
                    if not PWSCF_ecutwfc == "":
                        temp = float(re.sub("[^0-9.]", "", line).strip())
                        PWSCF_ecutrho = temp / PWSCF_ecutwfc
                    
                        print("PWSCF_ecutrho = "+str(PWSCF_ecutrho), file=summary)
                    else:
                        print("# WARN - No value for PWSCF_ecutwfc, PWSCF_ecutrho cannot be calculated for ["+refcode+"]", file=summary)
                        printToLog("# WARN - Compound ["+refcode+"] No value for PWSCF_ecutwfc, PWSCF_ecutrho cannot be calculated")

                elif "conv_thr" in line:
                    PWSCF_conv_thr = line[line.find("=")+1:].strip()

                    print("PWSCF_conv_thr = "+str(PWSCF_conv_thr), file=summary)            
                elif "CELL_PARAMETERS {alat}" in line:
                    temp = re.sub('\s{2,}', ' ', lines[number]).strip().split(" ")
                    reduction_factor = temp[0]

                    vector1 = np.array(re.sub('\s{2,}', ' ', lines[number]).strip().split(" "), dtype=float)
                    vector2 = np.array(re.sub('\s{2,}', ' ', lines[number+1]).strip().split(" "), dtype=float)
                    vector3 = np.array(re.sub('\s{2,}', ' ', lines[number+2]).strip().split(" "), dtype=float)

                    printToLog("# INFO - Compound ["+refcode+"] Reduced by a factor of ["+str(1/float(reduction_factor))+"]")
                elif "ATOMIC_POSITIONS" in line:  
                    start = number
                elif "K_POINTS automatic" in line:  
                    end = number - 3
            sub = lines[start:end]
            
            curr = 0
            symmetryEquivelents = []
            for number, line in enumerate(sub, 0):
                if number == curr:
                    equivalents = []
                    coordinates = re.sub('\s{2,}', ' ', line).strip().split(" ")
                    
                    x = float(coordinates[1])
                    y = float(coordinates[2])
                    z = float(coordinates[3])
    
                    for operation in CIF_symmetryElements:
                        operation = operation.split(",")
    
                        equivalent = []
                        equivalent.append(round(float(eval(str(operation[0]))), 5) % 1)
                        equivalent.append(round(float(eval(str(operation[1]))), 5) % 1)
                        equivalent.append(round(float(eval(str(operation[2]))), 5) % 1)      
                        
                        equivalents.append(str(equivalent))

                    unique = round(len(list(set(equivalents))) * float(reduction_factor))
                    symmetryEquivelents.append(unique)
                    curr += unique            
    else:
        printToLog("# WARN - Compound ["+refcode+"] does not have a .in file")
        quit()
            
    #REFCODE.out
    pwscfOut = os.path.join(refcodeDirectory, refcode+".out")
    if os.path.isfile(pwscfOut):
        with open(pwscfOut) as file:
            print("\n# -PWSCF output-\n", file=summary)
            
            lines = file.readlines()
            for number, line in enumerate(lines, 0):                  
                if "Program PWSCF" in line:
                    year = re.sub('\s{2,}', ' ', line.strip()).strip().split(" ")[5]
                    if len(year) < 9:
                        year = "0"+str(year)
                        
                    date = datetime.datetime.strptime(year+line.strip()[-8:].replace(" ", "0"), "%d%b%Y%H:%M:%S")
                    PWSCF_start_time = date.strftime("%Y-%m-%d %H:%M:%S")
                    PWSCF_version = line.strip().split(" ")[2]

                    print("PWSCF_version = "+str(PWSCF_version), file=summary)
                    print("PWSCF_start_time = "+str(PWSCF_start_time), file=summary)   
                elif "Number of MPI processes" in line:
                    PWSCF_numberMPI = float(re.sub("[^0-9.]", "", line).strip())

                    print("PWSCF_numberMPI = "+str(PWSCF_numberMPI), file=summary)
                elif "Threads/MPI process" in line:
                    PWSCF_numberThreads =  float(re.sub("[^0-9.]", "", line).strip())

                    print("PWSCF_numberThreads = "+str(PWSCF_numberThreads), file=summary)
                elif "R & G space division" in line:
                    PWSCF_RG = float(re.sub("[^0-9.]", "", line).strip())

                    print("PWSCF_RG = "+str(PWSCF_RG), file=summary)
                elif "Estimated total dynamical RAM" in line:
                    PWSCF_estimatedRAM = float(re.sub("[^0-9.]", "", line).strip())
                    
                    print("PWSCF_estimatedRAM = "+str(PWSCF_estimatedRAM), file=summary)
                elif "bfgs converged" in line:
                    PWSCF_scfCycles = float(re.sub("[^0-9.]", " ", line).strip()[:5].strip())
                    PWSCF_bfgsSteps = float(re.sub("[^0-9.]", " ", line).strip()[-5:].strip())
                    
                    print("PWSCF_scfCycles = "+str(PWSCF_scfCycles), file=summary)
                    print("PWSCF_bfgsSteps = "+str(PWSCF_bfgsSteps), file=summary) 
                elif "Final energy" in line:
                    PWSCF_finalEnergy = float(re.sub("[^0-9.-]", "", line).strip())
                    
                    print("PWSCF_finalEnergy = "+str(PWSCF_finalEnergy), file=summary)
                elif "This run was terminated on" in line:
                    temp = line.strip()[-9:].replace(" ", "0")+line.strip()[-19:-9].strip().replace(" ", "0")
                    date = datetime.datetime.strptime(temp, "%d%b%Y%H:%M:%S")
                    PWSCF_end_time = date.strftime("%Y-%m-%d %H:%M:%S")
                    
                    print("PWSCF_end_time = "+str(PWSCF_end_time), file=summary)      

                    temp = re.sub('\s{2,}', ' ', lines[number-3].strip()).strip()
                    if str(temp[0:5]) == "PWSCF":
                        PWSCF_CPU_time = temp.split(": ")[1].split(" CPU")[0]
                        PWSCF_WALL_time = temp.split("CPU ")[1].split(" WALL")[0]
                        print("PWSCF_CPU_time = "+str(PWSCF_CPU_time), file=summary)  
                        print("PWSCF_WALL_time = "+str(PWSCF_WALL_time), file=summary)
                    else:
                        printToLog("# WARN - No PWSCF CPU/WALL time found")   
                elif "JOB DONE" in line:
                    PWSCF_done = True
                    print("PWSCF_done = "+str(PWSCF_done), file=summary) 
        if PWSCF_done == False:
            print("PWSCF_done = "+str(PWSCF_done), file=summary) 
            print("# WARN - PWSCF did not run to completion", file=summary) 
            printToLog("# WARN - Compound ["+refcode+"] PWSCF did not run to completion")
            quit()

        if PWSCF_scfCycles == "":
            print("# WARN - Convergence not reached in PWSCF output", file=summary)
            printToLog("# WARN - Compound ["+refcode+"] Convergence not reached in PWSCF output")
        if PWSCF_finalEnergy == "":
            print("# WARN - Did not reach a final energy", file=summary)
            printToLog("# WARN - Compound ["+refcode+"] Did not reach a final energy")
    else:
        print("# WARN - No .out file found for compound with refcode ["+refcode+"]", file=summary)
        printToLog("# WARN - Compound ["+refcode+"] No PWSCF .out file found")
        quit()
        
    #gipaw.REFCODE.in
    gipawIn = os.path.join(refcodeDirectory, "gipaw."+refcode+".in")
    if os.path.isfile(gipawIn):
        print("\n# -GIPAW params-\n", file=summary)            
        with open(gipawIn) as file:
            lines = file.read().splitlines()
            for line in lines:
                if "q_gipaw" in line:
                    GIPAW_q_gipaw = float(re.sub("[^0-9.]", "", line).strip())

                    print("GIPAW_q_gipaw = "+str(GIPAW_q_gipaw), file=summary)
    else:
        print("# WARN - No gipaw .in file found for compound with refcode ["+refcode+"]", file=summary)
        printToLog("# WARN - Compound ["+refcode+"] No GIPAW .in file found")
        quit()

    #gipaw.REFCODE.out
    gipawOut = os.path.join(refcodeDirectory, "gipaw."+refcode+".out")
    if os.path.isfile(gipawOut):
        print("\n# -GIPAW output-\n", file=summary)            
        with open(gipawOut) as file:
            start = 0
            lines = file.readlines()
            for number, line in enumerate(lines, 0): 
                if "Program GIPAW" in line:
                    year = re.sub('\s{2,}', ' ', line.strip()).strip().split(" ")[5]
                    if len(year) < 9:
                        year = "0"+str(year)
                        
                    date = datetime.datetime.strptime(year+line.strip()[-8:].replace(" ", "0"), "%d%b%Y%H:%M:%S")
                    GIPAW_start_time = date.strftime("%Y-%m-%d %H:%M:%S")
                    GIPAW_version = line.strip().split(" ")[2]
                    
                    print("GIPAW_version = "+str(GIPAW_version), file=summary)
                    print("GIPAW_start_time = "+str(GIPAW_start_time), file=summary)   
                elif "Number of MPI processes" in line:
                    GIPAW_numberMPI = float(re.sub("[^0-9.]", "", line).strip())

                    print("GIPAW_numberMPI = "+str(GIPAW_numberMPI), file=summary)
                elif "Threads/MPI process" in line:
                    GIPAW_numberThreads = float(re.sub("[^0-9.]", "", line).strip())

                    print("GIPAW_numberThreads = "+str(GIPAW_numberThreads), file=summary)
                elif "R & G space division" in line:
                    GIPAW_RG = float(re.sub("[^0-9.]", "", line).strip())

                    print("GIPAW_RG = "+str(GIPAW_RG), file=summary)
                elif "Macroscopic shape contribution in ppm" in line:
                    GIPAW_mscPPM = float(re.sub("[^0-9.]", "", line).strip())

                    print("GIPAW_mscPPM = "+str(GIPAW_mscPPM), file=summary)
                elif "NMR macroscopic correction" in line:
                    GIPAW_msCorrection.append(re.sub('\s{2,}', ' ', lines[number+1]).strip().split(" "))
                    GIPAW_msCorrection.append(re.sub('\s{2,}', ' ', lines[number+2]).strip().split(" "))
                    GIPAW_msCorrection.append(re.sub('\s{2,}', ' ', lines[number+3]).strip().split(" "))

                    GIPAW_msCorrection = str(GIPAW_msCorrection)
                    print("GIPAW_msCorrection = "+str(GIPAW_msCorrection), file=summary)
                elif "Total sigma" in line:
                    if start == 0:
                        start = number
                        printToLog(f"# INFO - Compound [{refcode}] Sigma values start on line [{start}]")
                elif "This run was terminated on" in line:
                    temp = line.strip()[-9:].replace(" ", "0")+line.strip()[-19:-9].strip().replace(" ", "0")
                    date = datetime.datetime.strptime(temp, "%d%b%Y%H:%M:%S")
                    GIPAW_end_time = date.strftime("%Y-%m-%d %H:%M:%S")
                    
                    print("GIPAW_end_time = "+str(GIPAW_end_time), file=summary)

                    temp = re.sub('\s{2,}', ' ', lines[number-3].strip()).strip()
                    if str(temp[0:5]) == "GIPAW":
                        GIPAW_CPU_time = temp.split(": ")[1].split(" CPU")[0]
                        GIPAW_WALL_time = temp.split("CPU ")[1].split(" WALL")[0]
                        print("GIPAW_CPU_time = "+str(GIPAW_CPU_time), file=summary)  
                        print("GIPAW_WALL_time = "+str(GIPAW_WALL_time), file=summary)
                    else:
                        printToLog("# WARN - No GIPAW CPU/WALL time found")             
                elif "JOB DONE" in line:
                    GIPAW_done = True
                    print("GIPAW_done = "+str(GIPAW_done), file=summary)
            if GIPAW_done == False:
                print("GIPAW_done = "+str(GIPAW_done), file=summary) 
                print("# WARN - GIPAW did not run to completion", file=summary) 
                printToLog("# WARN - Compound ["+refcode+"] GIPAW did not run to completion")
                quit()

            print("\n# -Sigma values-\n", file=summary)            
            sub = lines[start:]
            printToLog(f"# INFO - Compound [{refcode}] Contains [{len(symmetryEquivelents)}] unique atoms")
            printToLog(f"# INFO - Compound [{refcode}] Number of symmetry equivalents of each atom [{symmetryEquivelents}]")

            count = 0
            previous = -10
            print("#BEGIN_ATOMIC_POSITIONS", file=summary)
            
            for number, line in enumerate(sub, 0):
                if count < len(symmetryEquivelents):
                    if number == previous + (10 * (symmetryEquivelents[count])):
                        printToLog("# INFO - Expecting ["+str(symmetryEquivelents[count])+"] symmetry equivalent atoms")
                        regex = re.compile('[^a-zA-Z ]')
    
                        previous = {}
                        for i in range(symmetryEquivelents[count]):
                            temp = number - (10 * i)
                            printToLog("# INFO - "+str(sub[temp]).strip())
                            curr = sub[temp]
    
                            atom = str(curr.strip()[:13].strip())
                            sigma = float(curr.strip()[-15:].lstrip())
                            
                            previous[atom] = sigma
    
                        currentSum = 0

                        for activeAtom, activeSigma in previous.items():
                            currentSum += activeSigma
                            for atom, sigma in previous.items():
                                diff = activeSigma - sigma
                                if diff < -tolerance or diff > tolerance:
                                    printToLog("# WARN - Compound ["+refcode+"] has symmetry equivalent atoms ["+str(activeAtom)+"] and ["+str(atom)+"] outside tolerance ["+str(diff)+"]")
                        matrix = []
                        matrix.append(re.sub('\s{2,}', ' ', sub[number+1]).strip().split(" "))
                        matrix.append(re.sub('\s{2,}', ' ', sub[number+2]).strip().split(" "))
                        matrix.append(re.sub('\s{2,}', ' ', sub[number+3]).strip().split(" "))

                        sigma = []
                        sigma.append(re.sub('\s{2,}', ' ', sub[number+6]).strip().split()[3])
                        sigma.append(re.sub('\s{2,}', ' ', sub[number+7]).strip().split()[3])
                        sigma.append(re.sub('\s{2,}', ' ', sub[number+8]).strip().split()[3])
                        sigma = sorted(sigma)
                        
                        line = re.sub('\s{2,}', ' ', line)
                        line = re.sub("Atom", "", line).lstrip().strip()
                        line = re.sub("pos: ", "", line)
                        line = re.sub("Total sigma: ", "", line)

                        coords = line.split("(")[1].split(")")[0]
                        line = line.split()
                        
                        print(f"#[{line[0]}] {line[1]} {line[-1]} ({str(round(currentSum / len(previous),2))}) SEP:{str(len(previous))} [XYZ{coords}] {matrix} [sigma_11 {sigma[0]}] [sigma_22 {sigma[1]}] [sigma_33 {sigma[2]}]",file=summary)
                        count += 1
                        previous = number
            print("#END_ATOMIC_POSITIONS", file=summary)
    else:
        print("# WARN - No gipaw .out file found for compound with refcode ["+refcode+"]", file=summary)
        printToLog("# WARN - Compound ["+refcode+"] No GIPAW .out file found")
        quit()

    printToLog("# INFO - Compound ["+refcode+"] Processing .cif file")
    with open(cif) as file:
        lines = file.readlines()
        
        cell_params = {}
        reduced_cell_params = {}
        print("\n# -CIF params-\n", file=summary)

        for number, line in enumerate(lines, 0):
            if "_cell_length_a" in line:
                cell_params["a"] = float(lines[number].split()[1].split("(")[0])
                cell_params["b"] = float(lines[number+1].split()[1].split("(")[0])
                cell_params["c"] = float(lines[number+2].split()[1].split("(")[0])
                cell_params["α"] = float(lines[number+3].split()[1].split("(")[0])
                cell_params["β"] = float(lines[number+4].split()[1].split("(")[0])
                cell_params["γ"] = float(lines[number+5].split()[1].split("(")[0])
                cell_params["volume"] = float(cellVolume(cell_params["a"], cell_params["b"], cell_params["c"], math.radians(cell_params["α"]), math.radians(cell_params["β"]), math.radians(cell_params["γ"])))
    
                for key, value in cell_params.items():
                    cell_params[key] = round(value,4)
                print("CIF_cell = ["+str(cell_params)+"]", file=summary)  
                printToLog("# INFO - Compound ["+refcode+"] Cell params ["+str(cell_params)+"]")

        temp = parseAlat(Alat, vector1, vector2, vector3)
        temp.append(float(cellVolume(float(temp[0]), float(temp[1]), float(temp[2]), math.radians(float(temp[3])), math.radians(float(temp[4])), math.radians(float(temp[5])))))
        temp = [str(round(float(elem), 4)) for elem in temp]
        
        reduced_cell_params["a"] = temp[0]
        reduced_cell_params["b"] = temp[1]
        reduced_cell_params["c"] = temp[2]
        reduced_cell_params["α"] = temp[3]
        reduced_cell_params["β"] = temp[4]
        reduced_cell_params["γ"] = temp[5]
        reduced_cell_params["volume"] = temp[6]
        print("CIF_reduced_cell = ["+str(reduced_cell_params)+"]", file=summary)  
        printToLog("# INFO - Compound ["+refcode+"] Reduced cell params ["+str(reduced_cell_params)+"]")

    printToLog("# INFO - Compound ["+refcode+"] Processing .out file")
    with open(out) as file:                
        lines = file.readlines()
        for number, line in enumerate(lines, 1):      
            if "Begin final coordinates" in line:
                start = number + 2                    
            elif "End final coordinates" in line:
                end = number - 1
        if len(lines[start:end]) == 0:
            printToLog("# WARN - Compound ["+refcode+"] Does not have final coordinates")
            quit()

    creator = mol2Creator(log=log,directory=refcodeDirectory,refcode=refcode,cell_params=reduced_cell_params,atom_positions=lines[start:end],df=df)
    creator.create()
    
    if os.path.isfile(gipawOut) and os.path.isfile(pwscfOut):
        if not PWSCF_numberMPI == GIPAW_numberMPI:
            printToLog("# WARN - Compound ["+refcode+"] Number of MPI do not match between PWSCF and GIPAW outputs")
        if not PWSCF_numberThreads == GIPAW_numberThreads:
            printToLog("# WARN - Compound ["+refcode+"] Number of threads does not match between PWSCF and GIPAW outputs")

shutil.copyfile(summaryPath, summaryCopyPath)
printToLog("# INFO - Compound ["+ refcode +"] Copied summary file ["+ summaryCopyPath + "]")

outputPath = os.path.join(output_files, refcode)
removeDirectory(outputPath, "# INFO - Compound ["+ refcode +"] Cleaning existing output path at")
printToLog("# INFO - Compound ["+ refcode +"] Copied output path ["+str(refcodeDirectory)+"] to ["+str(outputPath)+"]")
shutil.copytree(refcodeDirectory, outputPath)
