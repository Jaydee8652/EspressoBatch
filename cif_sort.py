#cif_sort - Jacob Duddridge - v2

# Filters unwanted .cif files based on data from the CSD (provided as a .csv) and the completeness of the structural data
# within the .cif and available .UPF files

# Cannot fully filter co-crystals, compounds with "structure-formula mismatch" likely contain an unreported co-crystal,
# and can be filtered manually with visual inspection in Mercury. Structures of this type are filtered to their own file
# for convenience

# All processes are reported to cif_sort.log for debugging

#Imports
import os
import csv
import re
import shutil
import sys
import pandas as pd
import datetime
import time
from utils.generic_utils import printToLog as pl, createDirectory as cd, removeDirectory as rd

#Params - can be modified
rCap = 10
volumeCap = 6000

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)
def removeDirectory(path, text):
    rd(log, path, text)
            
#Counts for stats
csvFilterCount = 0
structureFormulaMismatchCount = 0
noStructureDataCount = 0
noHydrogenDataCount = 0
incompleteHydrogenDataCount = 0
unaccountedElementsCount = 0

unaccountedElements = {}

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

#Make sure there is a directory to sort
pathToSort = os.path.join(homeDirectory, "Original_CIFs")
createDirectory(pathToSort, "# WARN - No directory found for .cifs to sort. Place .cif files in or replace the newly created directory at", True)

cifFiles = [file for file in os.listdir(pathToSort) if file.endswith('.cif') and os.path.isfile(os.path.join(pathToSort, file))]#Get .cifs from directory
numberOfCifs = len(cifFiles)
survivingStructures = numberOfCifs

if numberOfCifs == 0:#Make sure there are .cifs in the directory
    printToLog("# WARN - No .cif files found to sort. Place .cif files in ["+ pathToSort + "]")
    sys.exit()
else:
    printToLog("# INFO - " + str(numberOfCifs) + " .cif files found at ["+ pathToSort + "]")

cifDict = {}
for file in cifFiles:#Put cifs into a dict for easy access
    fileName = os.path.splitext(file)[0]
    fileName = fileName.replace(".cif", "")#Cull extraneous file name data. Only necessary if the files were saved from CSD weirdly
    cifDict[fileName] = file

#Make sure there is a directory for PSEUDOS
pseudosPath = os.path.join(homeDirectory, "PSEUDOS/")
createDirectory(pseudosPath, "# WARN - No directory found for PSEUDOS. Place .UPF files in or replace the newly created directory at", True)

pseudosFiles = [file for file in os.listdir(pseudosPath) if file.endswith('.UPF') and os.path.isfile(os.path.join(pseudosPath, file))]#Get .UPFs from directory
numberOfPsueds = len(pseudosFiles) # determine number of .UPF files
psuedoElements = list()
if numberOfPsueds == 0:#Make sure there are .UPFs in the directory
    printToLog("# WARN - No .UPF files found. Place .UPF files in ["+ pseudosPath + "]")
    sys.exit()
else:
    printToLog("# INFO - " + str(numberOfPsueds) + " .UPF files found at ["+ pseudosPath + "]")

    psuedoElements = [file.split(".")[0] for file in pseudosFiles]
    printToLog("# INFO - Following elements accounted for ["+str(list(set(psuedoElements)))+"]")

#Make sure there are directories to put the sorted files
cifPath = os.path.join(homeDirectory, "cifs")
createDirectory(cifPath, "# INFO - No file found storing .cifs, created at", False)

#Make sure there are directories to put the sorted files
validatedPath = os.path.join(cifPath, "Validated_CIFs")
createDirectory(validatedPath, "# INFO - No file found for validated .cifs, created at", False)

discardedPath = os.path.join(cifPath, "Discarded_CIFs")
createDirectory(discardedPath, "# INFO - No file found for discarded .cifs, created at", False)

structureFormulaMismatchPath = os.path.join(discardedPath, "StructureFormulaMismatched_CIFs")
createDirectory(structureFormulaMismatchPath, "# INFO - No file found for structure-formula mismatched .cifs, created at", False)

noStructureDataPath = os.path.join(discardedPath, "NoStructureData_CIFs")
createDirectory(noStructureDataPath, "# INFO - No file found for empty .cifs, created at", False)

noHydrogenDataPath = os.path.join(discardedPath, "NoHydrogenData_CIFs")
createDirectory(noHydrogenDataPath, "# INFO - No file found for .cifs with no hydrogen data, created at", False)

incompleteHydrogenDataPath = os.path.join(discardedPath, "IncompleteHydrogenData_CIFs")
createDirectory(incompleteHydrogenDataPath, "# INFO - No file found for .cifs with incomplete hydrogen data, created at", False)

unaccountedElementsPath = os.path.join(discardedPath, "UnaccountedElements_CIFs")
createDirectory(unaccountedElementsPath, "# INFO - No file found for .cifs with unaccounted element types, created at", False)

#Make sure there is a .csv to read from
structureDataPath = os.path.join(homeDirectory, "structure_data.csv")
if not os.path.exists(structureDataPath):
    printToLog("# WARN - No .csv file found to load compound data. Copy .csv from the CSD to the following path ["+ structureDataPath + "]")
    sys.exit()

with open("structure_data.csv","r", encoding="utf-8-sig") as fileCSV:#Open accompanying .csv - Obtained from CSD
    readCSV = csv.DictReader(fileCSV)
    for line in readCSV:
        refcode = str(line["[REFCODE]"])
        cellVolume = float(line["[_cell_volume]"])
        rcellVolume = float(line["[_rcell_volume]"])
        rFactor = float(line["[_refine_ls_R_factor]"])
        disorder = line["[_exptl_DISORDER]"]

        printToLog("# INFO - Processing compound with refcode ["+ refcode +"]")
        if not cifDict.__contains__(refcode):#Warn if there is data but no .cif
            printToLog("# WARN - No .cif file found for structure with refcode [" + refcode + "]")
        else:
            fileName = cifDict.get(refcode)
            with open(os.path.join(pathToSort, fileName), "r") as cif:
                printToLog("# INFO - Processing .cif file for compound with refcode [" + refcode + "]")
                formulaDict = {}
                workingDict = {}
                expl = ""
                
                for line in cif:
                    line = line.strip()
                    if len(line) == 0 or line.startswith("#"):
                        if "#Symmetry data" in line:                        
                            appended = True
                        else:
                            continue#   skip blank lines and comments
                    if m := re.match(r"_?(\w+)\s+(.*)", line):
                        key = m.group(1)
                        value = m.group(2)
                        
                        if key.startswith("chemical_formula_sum"):
                            formula = value.strip("\'")
                            printToLog("# INFO - Compound [" + refcode + "] found to have formula ["+ formula +"]")
                            regex = re.compile('[^a-zA-Z ]')
                            atoms = regex.sub('', formula).strip().split(" ")
                            atoms = list(set([atom.lower().capitalize() for atom in atoms]))

                            unaccounted = []
                            printToLog("# INFO - Compound [" + refcode + "] found to contain the following atom types: ["+ str(atoms) +"]")
                            for atom in atoms:
                                if not psuedoElements.__contains__(atom):
                                    unaccounted.append(atom)
                                    #printToLog("# WARN - Compound [" + refcode + "] found to contain atom not accounted for by .UPF files: ["+ str(atom) +"]")
                                    if atom in unaccountedElements:
                                        unaccountedElements[atom] += 1
                                    else:
                                        unaccountedElements[atom] = 1
                            formula = formula.split()
                            for element in formula:#Save formula to a dict
                                formulaDict[re.findall('\\d+|\\D+', element)[0]] = int(re.findall('\\d+|\\D+', element)[1])
                            workingDict = formulaDict.copy()
                        else:
                            elementSymbolLine = re.findall('\\d+|\\D+', value)[0].strip(" -")
                            if "?" in value and expl == "":
                                expl = expl + "Disordered [" + str(value) + "]"
                            elif workingDict.__contains__(elementSymbolLine):
                                workingDict[elementSymbolLine] -= 1#Increment a copy of the formula dict for every corresponding atom position, counts down for convenience
                                
                #Generate explanation for discard
                if cellVolume > volumeCap:
                    expl = expl + "Cell Volume [" + str(cellVolume) + "] greater than volume cap [" + str(volumeCap) + "]"
                if rFactor > rCap:
                    if not expl == "":
                        expl += " & "
                    expl = expl + "R factor [" + str(rFactor) + "] greater than cutoff point [" + str(rCap) + "]"
                if not disorder == "":
                    if not expl == "":
                        expl += " & "
                    expl = expl + "Disordered [" + str(disorder) + "]"
    
                #Filter by volume, r factor and disorder
                if not expl == "":
                    shutil.copyfile(os.path.join(pathToSort, fileName), os.path.join(discardedPath, refcode + ".cif"))
                    printToLog("# INFO - Compound [" + refcode + "] discarded. " + expl)
                    csvFilterCount += 1
                    survivingStructures -= 1
                    continue
                if not all(workingDict.get(key) % formulaDict.get(key) == 0 and not workingDict.get(key) == formulaDict.get(key) for key in formulaDict):#If the formula doesn't equal 0 (All atoms accounted for), or a negative multiple (Symmetry equivalent), discard TODO: Cannot account for co-crytsals.
                    if all(workingDict.get(key) == formulaDict.get(key) for key in formulaDict):
                        printToLog("# INFO - Compound [" + refcode + "] discarded. No structural data, all atoms unaccounted for: ["+ str(workingDict) +"]")
                        shutil.copyfile(os.path.join(pathToSort, fileName), os.path.join(noStructureDataPath, refcode + ".cif"))

                        noStructureDataCount += 1
                    elif (workingDict.get("H") == formulaDict.get("H")):
                        printToLog("# INFO - Compound [" + refcode + "] discarded. No hydrogen data, unaccounted for atoms: ["+ str(workingDict) +"]")
                        shutil.copyfile(os.path.join(pathToSort, fileName), os.path.join(noHydrogenDataPath, refcode + ".cif"))

                        noHydrogenDataCount += 1
                    elif all(workingDict.get(key) == 0 or key == "H" for key in formulaDict):
                        printToLog("# INFO - Compound [" + refcode + "] discarded. Incomplete hydrogen data, unaccounted for atoms: ["+ str(workingDict) +"]")
                        shutil.copyfile(os.path.join(pathToSort, fileName), os.path.join(incompleteHydrogenDataPath, refcode + ".cif"))

                        incompleteHydrogenDataCount += 1
                    else:
                        printToLog("# INFO - Compound [" + refcode + "] discarded. Structure-formula mismatch, unaccounted for atoms: ["+ str(workingDict) +"]")
                        shutil.copyfile(os.path.join(pathToSort, fileName), os.path.join(structureFormulaMismatchPath, refcode + ".cif"))
                        
                        structureFormulaMismatchCount += 1
                    survivingStructures -= 1
                else:
                    if len(unaccounted) > 0:
                        printToLog("# INFO - Compound [" + refcode + "] discarded. Found to contain the following atoms not accounted for by .UPF files: ["+ str(unaccounted) +"]")
                        shutil.copyfile(os.path.join(pathToSort, fileName), os.path.join(unaccountedElementsPath, refcode + ".cif"))

                        unaccountedElementsCount += 1
                        survivingStructures -= 1
                    else:
                        printToLog("# INFO - Compound [" + refcode + "] validated")
                        shutil.copyfile(os.path.join(pathToSort, fileName), os.path.join(validatedPath, refcode + ".cif"))



finalPath = os.path.join(cifPath, "Original_CIFs")
removeDirectory(finalPath, "# INFO - Cleaning existing sorted path at")
printToLog("# INFO - Moving sorted path ["+str(pathToSort)+"] to ["+str(finalPath)+"]")
shutil.move(pathToSort, finalPath)

printToLog("# INFO - [" + str(csvFilterCount) + "] compounds filtered by .csv properties")
printToLog("# INFO - [" + str(noStructureDataCount) + "] compounds found to have no structural data")
printToLog("# INFO - [" + str(noHydrogenDataCount) + "] compounds found to have no hydrogen data")
printToLog("# INFO - [" + str(incompleteHydrogenDataCount) + "] compounds found to have incomplete hydrogen data")
printToLog("# INFO - [" + str(structureFormulaMismatchCount) + "] compounds found to have structure-formula mismatch")
printToLog("# INFO - [" + str(unaccountedElementsCount) + "] compounds found to have atoms unaccounted for by .UPF files")

if len(unaccountedElements) > 0:
    printToLog("# WARN - The following atoms are unaccounted for by the current PSUEDOS: ["+str(unaccountedElements)+"]")
else:
    printToLog("# INFO - All present atoms are accounted for by the current PSEUDOS")
survivalPercent = round((survivingStructures / numberOfCifs) * 100, 3)
printToLog("# INFO - Sort complete, [" + str(survivingStructures) + "] of [" + str(numberOfCifs) + "] compounds validated, with a % survival of [" +str(survivalPercent)+ "%]")
