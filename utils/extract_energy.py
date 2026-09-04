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

#Jank thing to fix the path. Very annoying artefact of running python scripts by absolute path.
sys.path[0] = sys.path[0][:-6] + sys.path[0][-6:].replace("/utils", "")

from utils.generic_utils import printToLog as pl, createDirectory as cd, removeDirectory as rd, writeCSV

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)    
def createDirectory(path, text, exit):
    cd(log, path, text, exit)
def removeDirectory(path, text, exit):
    rd(log, path, text)

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")

refcodeDirectory = os.getcwd() #Directory where we are
input_path = os.path.split(refcodeDirectory)[0]
homeDirectory = os.path.split(input_path)[0]

refcode = os.path.basename(refcodeDirectory)
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Compound ["+refcode+"] Starting "+str(os.path.basename(sys.argv[0]).split(".")[0])+" in ["+str(refcodeDirectory)+"]")

# Get .cif
cifPath = os.path.join(refcodeDirectory, refcode+".cif")
if os.path.isfile(cifPath):
    with open(cifPath, "r") as cif:
        for line in cif: 
            if "_cell_formula_units_Z" in line:
                cell_formula_units_Z = float(re.sub("[^0-9.]", "", line).strip())
            elif "_symmetry_space_group_name" in line:
                symmetry_space_group_name = str(' '.join(line.split()[1:])).strip("'")
else:
    printToLog("# WARN - Compound ["+refcode+"] does not have a .cif file")
    quit()

# Remove .save directory
saveDirectory = os.path.join(refcodeDirectory, refcode+".save")
if os.path.exists(saveDirectory):
    printToLog("# INFO - Compound ["+ refcode +"] Removing .save file ["+ saveDirectory + "]")
    shutil.rmtree(saveDirectory)

# Remove .wfc files
printToLog("# INFO - Compound ["+ refcode +"] Removing .wfc files")
regex = re.compile('[^a-zA-Z.]')
wfcFiles = [file for file in os.listdir(refcodeDirectory) if regex.sub('', file).endswith('.wfc') and os.path.isfile(os.path.join(refcodeDirectory, file))]
for wfc in wfcFiles:
    os.remove(os.path.join(refcodeDirectory, wfc))

# Remove .mix files
printToLog("# INFO - Compound ["+ refcode +"] Removing .mix files")
regex = re.compile('[^a-zA-Z.]')
mixFiles = [file for file in os.listdir(refcodeDirectory) if regex.sub('', file).endswith('.mix') and os.path.isfile(os.path.join(refcodeDirectory, file))]
for mix in mixFiles:
    os.remove(os.path.join(refcodeDirectory, mix))

# Clean existing summary
summaryPath = os.path.join(refcodeDirectory, refcode+"_summary.txt")
if os.path.exists(summaryPath):
    printToLog("# INFO - Compound ["+ refcode +"] Cleaning existing summary file ["+ summaryPath + "]")
    os.remove(summaryPath)
    
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
            if not read.__contains__("BATCH_done"):
                print("BATCH_end_time = "+str(now), file=file)       
                print("BATCH_done = "+str(True), file=file)            
        with open(batch) as file:
            lines = file.read().splitlines()
            for line in lines:
                print(line, file=summary)
    else:
        print("# WARN - No batch file found for compound with refcode ["+refcode+"]", file=summary)
        printToLog("# WARN - Compound ["+refcode+"] No batch file found")

    #REFCODE.out
    out = os.path.join(refcodeDirectory, refcode+".out")
    if os.path.isfile(out):
        with open(out) as file:
            print("\n# -scf output-\n", file=summary)
            
            lines = file.readlines()
            reduction_factor = 1
            for number, line in enumerate(lines, 1):
                if "crystal axes: (cart. coord. in units of alat)" in line:
                    temp = lines[number].split("=")[1].split("(")[1].split(")")[0]
                    temp = re.sub('\s{2,}', ' ', temp).strip().split(" ")
                    
                    reduction_factor = temp[0]
                    printToLog("# INFO - Compound ["+refcode+"] Reduced by a factor of ["+str(1/float(reduction_factor))+"]")
                if line.startswith("!"):
                    SCF_final_energy_ry = float(line[line.find("=")+1:].lstrip().split()[0]) * (1 / float(reduction_factor))
                    SCF_final_energy_kJ_mol1_molecule1 = 6.02214076e+23 * 2.1798741e-21 * (1 / int(cell_formula_units_Z)) * float(SCF_final_energy_ry)

                    print("CIF_cell_formula_units_Z = "+str(cell_formula_units_Z), file=summary)  
                    print("CIF_symmetry_space_group_name = "+str(symmetry_space_group_name), file=summary)  
                    print("CIF_reduction_factor = "+str(1 / float(reduction_factor)))  
                    print("SCF_final_energy_ry = "+str(SCF_final_energy_ry), file=summary)  
                    print("SCF_final_energy_kJ_mol⁻¹_molecule⁻¹ = "+str(SCF_final_energy_kJ_mol1_molecule1), file=summary)  
    else:
        print("# WARN - No .out file found for compound with refcode ["+refcode+"]", file=summary)
        printToLog("# WARN - Compound ["+refcode+"] No PWSCF .out file found")

localSheet = os.path.join(input_path, "sanity_sheet.csv")
if not os.path.isfile(localSheet):
    printToLog("# WARN - No local .csv found. Attempting to initialise")
    with open(localSheet, 'a') as file:
        file.write("[REFCODE]")

df = pd.read_csv(localSheet)  
df.set_index('[REFCODE]', inplace = True)
df = df.astype(object)

with open(summaryPath, "r") as file:            
    read = file.read()
    lines = read.splitlines()
    
    for line in lines:
        if not len(line) == 0 and not line.startswith("#"):
            value = line[line.find("=")+1:].strip()
            name = line[:line.find("=")-1].strip()
            writeCSV(df, refcode, "["+str(name)+"]", str(value))

df = df.replace("nan", "")
df.to_csv(localSheet) # Update local csv
printToLog("# INFO - Compound ["+refcode+"] Appended to .csv")
