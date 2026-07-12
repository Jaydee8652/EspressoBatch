# post_processing.py - Jacob Duddridge

# Extracts desired data from PWSCF and gipaw .out files to a summary file, saved here and in a dedicated directory
# Automatically removes symmetry equivalent atoms from the gipaw output

# Automatically run by slurm within the input directory

# All processes are reported to post_processing.log for debugging, all post_processing scripts output to the same log within a home directory,
# and individual processes produce their own log locally within the input file

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
from wyckoff import WyckoffDatabase

#Params - can be modified
tolerance = 0.01

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

#Create directory if it doesn't exist. Optionally crash deliberately if doesn't exist
def createDirectory(path, text, exit):
    if not os.path.exists(path):
        printToLog(text + " ["+ path + "]")
        os.makedirs(path)
        if exit:
            quit()

#Remove directory is it exists
def removeDirectory(path, text):
    if os.path.exists(path):
        printToLog(text + " ["+ path + "]")
        shutil.rmtree(path)

CIF_symmetryTableNumber = ""
CIF_symmetryElements = []
CIF_symmetryFactor = ""

PWSCF_ecutwfc = ""
PWSCF_ecutrho = "" 
PWSCF_conv_thr = ""
PWSCF_version = ""
PWSCF_time = ""
PWSCF_numberMPI = ""
PWSCF_numberThreads = ""
PWSCF_RG = ""
PWSCF_estimatedRAM = ""
PWSCF_scfCycles = ""
PWSCF_bfgsSteps = ""
PWSCF_finalEnergy = ""

PWSCF_done = False

GIPAW_q_gipaw = ""
GIPAW_version = ""
GIPAW_time = ""
GIPAW_numberMPI = ""
GIPAW_numberThreads = ""
GIPAW_RG = ""
GIPAW_mscPPM = ""
GIPAW_msCorrection = []

GIPAW_done = False


#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")

#homeDirectory = sys.argv[1]
refcodeDirectory = os.getcwd()#Directory where we are
homeDirectory = os.path.split(os.path.split(refcodeDirectory)[0])[0]

refcode = os.path.basename(refcodeDirectory)
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Compound ["+refcode+"] Starting "+str(os.path.basename(sys.argv[0]).split(".")[0])+" in ["+str(refcodeDirectory)+"]")

summariesPath = os.path.join(homeDirectory, "Summary_Files")
createDirectory(summariesPath, "# INFO - No directory found for summary files, created at ", False)

outputsPath = os.path.join(homeDirectory, "Output_Files")
createDirectory(outputsPath, "# INFO - No directory found for output files, created at", False)

# Get .cif directory
cifPath = os.path.join(refcodeDirectory, refcode+".cif")
if os.path.isfile(cifPath):
    with open(cifPath, "r") as cif:
        for line in cif:
            if line[0].isdigit():
                CIF_symmetryElements.append(line.split()[1])
            elif "_symmetry_Int_Tables_number" in line:
                CIF_symmetryTableNumber = int(re.sub("[^0-9.]", "", line).strip())         
else:
    printToLog("# WARN - Compound ["+refcode+"] does not have a .cif file")
    quit()

WyckoffData = WyckoffDatabase().data[CIF_symmetryTableNumber]

# Get structure data .csv
structureDataPath = os.path.join(homeDirectory, "structure_data.csv")
if os.path.exists(structureDataPath):
    df = pd.read_csv(structureDataPath, encoding="utf-8-sig")
    df.set_index('[REFCODE]', inplace = True)
    cellVolume = float(df.at[refcode, "[_cell_volume]"])
    rcellVolume = float(df.at[refcode, "[_rcell_volume]"])
    CIF_symmetryFactor = round(len(CIF_symmetryElements) / round(cellVolume / rcellVolume))
else:
    printToLog("# WARN - No .csv file found to load compound data. Copy .csv from the CSD to the following path ["+ structureDataPath + "]")
    quit()

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
        with open(batch) as file:
            lines = file.read().splitlines()
            for line in lines:
                print(line, file=summary)
    else:
        print("# WARN - No batch file found for compound with refcode ["+refcode+"]", file=summary)
        printToLog("# WARN - Compound ["+refcode+"] No batch file found")

    print("\n# -CIF data\n", file=summary)
    print("CIF_symmetryTableNumber = "+str(CIF_symmetryTableNumber), file=summary)
    print("CIF_symmetryElements = "+str(CIF_symmetryElements), file=summary)
    print("CIF_symmetryFactor = "+str(CIF_symmetryFactor), file=summary)
    
    #REFCODE.in
    pwscfIn = os.path.join(refcodeDirectory, refcode+".in")
    if os.path.isfile(pwscfIn):
        print("\n# -PWSCF params-\n", file=summary)            
        with open(pwscfIn) as file:
            lines = file.read().splitlines()
            for line in lines:
                if "ecutwfc" in line:
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
    else:
        print("# WARN - No .in file found for compound with refcode ["+refcode+"]", file=summary)
        printToLog("# WARN - Compound ["+refcode+"] No .in file found")
    
    #REFCODE.out
    pwscfOut = os.path.join(refcodeDirectory, refcode+".out")
    if os.path.isfile(pwscfOut):
        with open(pwscfOut) as file:
            print("\n# -PWSCF output-\n", file=summary)

            lineNumber = 0
            read = file.read()
            lines = read.splitlines()

            for line in lines:
                lineNumber += 1
                if "Program PWSCF" in line:
                    year = re.sub('\s{2,}', ' ', line.strip()).strip().split(" ")[5]
                    if len(year) < 9:
                        year = "0"+str(year)
                        
                    date = datetime.datetime.strptime(year+line.strip()[-8:].replace(" ", "0"), "%d%b%Y%H:%M:%S")
                    PWSCF_time = date.strftime("%Y-%m-%d %H:%M:%S")
                    PWSCF_version = line.strip().split(" ")[2]

                    print("PWSCF_version = "+str(PWSCF_version), file=summary)
                    print("PWSCF_time = "+str(PWSCF_time), file=summary)   
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
                elif "JOB DONE" in line:
                    PWSCF_done = True
                    print("PWSCF_done = "+str(PWSCF_done), file=summary) 
        if PWSCF_done == False:
            print("PWSCF_done = "+str(PWSCF_done), file=summary) 
            print("# WARN - PWSCF did not run to completion", file=summary) 

            printToLog("# WARN - Compound ["+refcode+"] PWSCF did not run to completion")
        if PWSCF_scfCycles == "":
            printToLog("# WARN - Compound ["+refcode+"] Convergence not reached in PWSCF output")
        if PWSCF_finalEnergy == "":
            printToLog("# WARN - Compound ["+refcode+"] Did not reach a final energy")
    else:
        print("# WARN - No .out file found for compound with refcode ["+refcode+"]", file=summary)
        printToLog("# WARN - Compound ["+refcode+"] No PWSCF .out file found")

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
        
    #gipaw.REFCODE.out
    gipawOut = os.path.join(refcodeDirectory, "gipaw."+refcode+".out")
    sigmaStart = 0
    if os.path.isfile(gipawOut):
        print("\n# -GIPAW output-\n", file=summary)            
        with open(gipawOut) as file:
            lineNumber = 0
            lines = file.read().splitlines()
            for line in lines:
                lineNumber += 1
                
                if "Program GIPAW" in line:
                    year = re.sub('\s{2,}', ' ', line.strip()).strip().split(" ")[5]
                    if len(year) < 9:
                        year = "0"+str(year)
                        
                    date = datetime.datetime.strptime(year+line.strip()[-8:].replace(" ", "0"), "%d%b%Y%H:%M:%S")
                    GIPAW_time = date.strftime("%Y-%m-%d %H:%M:%S")
                    GIPAW_version = line.strip().split(" ")[2]
                    
                    print("GIPAW_version = "+str(GIPAW_version), file=summary)
                    print("GIPAW_time = "+str(GIPAW_time), file=summary)   
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
                    GIPAW_msCorrection.append(re.sub('\s{2,}', ' ', lines[lineNumber]).strip().split(" "))
                    GIPAW_msCorrection.append(re.sub('\s{2,}', ' ', lines[lineNumber+1]).strip().split(" "))
                    GIPAW_msCorrection.append(re.sub('\s{2,}', ' ', lines[lineNumber+2]).strip().split(" "))

                    GIPAW_msCorrection = str(GIPAW_msCorrection)
                    print("GIPAW_msCorrection = "+str(GIPAW_msCorrection), file=summary)
                elif "Total sigma" in line:
                    if sigmaStart == 0:
                        sigmaStart = lineNumber                    
                    if ((lineNumber - sigmaStart) / 10) % CIF_symmetryFactor == 0:
                        regex = re.compile('[^a-zA-Z ]')
                        initialAtom = str(line.strip()[:13].strip())
                        initialElement = regex.sub('', initialAtom).strip().split(" ")[-1]
                        
                        previousDict = {}                        
                        for i in range(CIF_symmetryFactor):
                            curr = lines[(10 * i) + lineNumber - 1]

                            atom = str(curr.strip()[:13].strip())
                            sigma = float(curr.strip()[-15:].lstrip())
                            previousDict[atom] = sigma

                            coordinates = re.sub('\s{2,}', ' ', curr[curr.find("(")+1:curr.find(")")]).strip().split(" ")
                            x = float(coordinates[0])
                            y = float(coordinates[1])
                            z = float(coordinates[2])

                            for wyckoff in WyckoffData.wyckoff_positions:
                                #if not wyckoff.multiplicity == int(CIF_symmetryFactor):
                                for position in wyckoff.coordinates:                                    
                                    target = []
                                    target.append(float(eval(str(position[0]))))
                                    target.append(float(eval(str(position[1]))))
                                    target.append(float(eval(str(position[2]))))

                                    printToLog(atom +" "+ str(coordinates) + " " + str(target) + " "+ str(position))
                                    if all(round(float(target[i]) - float(coordinates[i]), 5) == 0 for i in range(len(target))):
                                    #        printToLog(atom +" "+ str(coordinates) + " " + str(target) + " "+ str(position))
                                            printToLog("HI!")
                                                #sigmaStart -= 10 * (CIF_symmetryFactor / wyckoff.multiplicity)
                                     #   for i in range(len(target)):
                                     #       printToLog(round(float(target[i]) - float(coordinates[i]), 5)))
                                            
                        #        for symmetryElement in CIF_symmetryElements:
                    #            if symmetryElement == "x,y,z":
                     #               continue
                                    
                      #          functions = symmetryElement.split(",")
                       #         printToLog(str(functions))
                                
                       #         x = float(position[0])
                       #         y = float(position[1])
                        #        z = float(position[2])

                         #       newPosition = []
                          #      newPosition.append(eval(functions[0]))
                          #      newPosition.append(eval(functions[1]))
                          #      newPosition.append(eval(functions[2]))

                        #        for i in range(len(position)):
                         #           printToLog(str(round(float(position[i]), 7) - round(float(newPosition[i]))))
                                
                        #        if all(round(float(position[i]), 7) - round(float(newPosition[i]), 7) == 0 for i in range(len(position))):
                           #         printToLog("HI!--------------------------------------------------------------------------------------------------------------------")#

                            
                            
                      #      if position == ['0.000000', '0.000000', '0.000000']:
                      #          printToLog("# WARN - Compound ["+refcode+"] Atom on inversion centre ["+str(atom)+"] ["+str(position)+"]")
                       #         previousDict.clear()
                       #         break
#
                 #       if len(previousDict) == 0:
                #            for i in range(round(CIF_symmetryFactor / 2)):
                #                curr = lines[(10 * i) + lineNumber - 1]
    
                #                atom = str(curr.strip()[:13].strip())
                 #               sigma = float(curr.strip()[-15:].lstrip())
                  #              previousDict[atom] = sigma
                   #             sigmaStart -= 10
                        
                        workingDict = previousDict.copy()
                        for atom in previousDict.keys():
                            element = regex.sub('', atom).strip().split(" ")[-1]
                            
                            if not element == initialElement and atom in workingDict:
                                del workingDict[atom]
                                sigmaStart -= 10
                                printToLog("# WARN - Compound ["+refcode+"] Symmetry equivalent atom of different type ["+str(initialAtom)+"] and ["+str(atom)+"]")

                        currentSum = 0
                        for activeAtom, activeSigma in workingDict.items():
                            currentSum += activeSigma
                            for atom, sigma in workingDict.items():
                                diff = activeSigma - sigma
                                if diff < -tolerance or diff > tolerance:
                                    printToLog("# WARN - Compound ["+refcode+"] has symmetry equivalent atoms ["+str(activeAtom)+"] and ["+str(atom)+"] outside tolerance ["+str(diff)+"]")
                        array = []
                        array.append(re.sub('\s{2,}', ' ', lines[lineNumber]).strip().split(" "))
                        array.append(re.sub('\s{2,}', ' ', lines[lineNumber+1]).strip().split(" "))
                        array.append(re.sub('\s{2,}', ' ', lines[lineNumber+2]).strip().split(" "))
                        array = str(array)
                        
                        print("_"+ str(line.lstrip()) + " - Averaged over "+str(len(workingDict))+" atoms: "+str(currentSum / len(workingDict)), file=summary)
                        print("_"+array, file=summary)

                        
                elif "JOB DONE" in line:
                    GIPAW_done = True
                    print("GIPAW_done = "+str(GIPAW_done), file=summary)
            if GIPAW_done == False:
                print("GIPAW_done = "+str(GIPAW_done), file=summary) 
                print("# WARN - GIPAW did not run to completion", file=summary) 
                printToLog("# WARN - Compound ["+refcode+"] GIPAW did not run to completion")
    else:
        print("# WARN - No gipaw .out file found for compound with refcode ["+refcode+"]", file=summary)
        printToLog("# WARN - Compound ["+refcode+"] No GIPAW .out file found")

    if os.path.isfile(gipawOut) and os.path.isfile(pwscfOut):
        if not PWSCF_numberMPI == GIPAW_numberMPI:
            printToLog("# WARN - Compound ["+refcode+"] Number of MPI do not match between PWSCF and GIPAW outputs")
        if not PWSCF_numberThreads == GIPAW_numberThreads:
            printToLog("# WARN - Compound ["+refcode+"] Number of threads does not match between PWSCF and GIPAW outputs")
            
summaryCopyPath =os.path.join(summariesPath, refcode+"_summary.txt")
if os.path.exists(summaryCopyPath):
    printToLog("# INFO - Compound ["+ refcode +"] Cleaning copied summary file ["+ summaryCopyPath + "]")
    os.remove(summaryCopyPath)
shutil.copyfile(summaryPath, summaryCopyPath)
printToLog("# INFO - Compound ["+ refcode +"] Copied summary file ["+ summaryCopyPath + "]")

outputPath = os.path.join(outputsPath, refcode)
removeDirectory(outputPath, "# INFO - Compound ["+ refcode +"] Cleaning existing output path at")
printToLog("# INFO - Compound ["+ refcode +"] Copied output path ["+str(refcodeDirectory)+"] to ["+str(outputPath)+"]")
shutil.copytree(refcodeDirectory, outputPath)
