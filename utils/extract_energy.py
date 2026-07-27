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
from generic_utils import writeCSV

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    logs = os.path.join(homeDirectory, "logs")
    if not os.path.exists(logs):
        os.makedirs(logs)
    info = str(info)

    time = ""
    if not info.startswith(" ---"):
        time = str(datetime.datetime.now().strftime("[%H:%M:%S] "))
    print(time+str(info))
    with open(log, "a") as l:
        l.write(time+str(info) + "\n")        
    with open(os.path.join(logs, log), "a") as l:
        l.write(time+str(info) + "\n")

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
else:
    printToLog("# WARN - Compound ["+refcode+"] does not have a .cif file")
    quit()

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
            for number, line in enumerate(lines, 1):                  
                if line.startswith("!"):
                    SCF_final_energy_ry = line[line.find("=")+1:].lstrip().split()[0]
                    SCF_final_energy_kcal_mol1_molecule1 = 5.2065398394955 * (10 ** -22) * 6.02214076 * (10 ** +23) * (1 / int(cell_formula_units_Z)) * float(SCF_final_energy_ry)

                    print("SCF_final_energy_ry = "+str(SCF_final_energy_ry), file=summary)  
                    print("SCF_final_energy_kcal_mol1_molecule1 = "+str(SCF_final_energy_kcal_mol1_molecule1), file=summary)  
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
df.to_csv(localSheet)#Update local csv
printToLog("# INFO - Compound ["+refcode+"] Appended to .csv")
