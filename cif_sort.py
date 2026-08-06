#cif_sort - Jacob Duddridge - v3

# Filters unwanted .cif files based on data from the CSD (provided as a .csv) and the completeness of the structural data
# within the .cif and available .UPF files

# All processes are reported to cif_sort.log for debugging

#Imports
import os
import re 
import pandas as pd
import io
import csv
import datetime
import time       
import sys
import shutil
import math
import numpy as np
from utils.generic_utils import printToLog as pl, createDirectory as cd, removeDirectory as rd, cellVolume
from utils.git_utils import downloadCSV, uploadCSV, appendCSV, updateCSV, batchCalculations, verify

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
            
#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

#Make sure there is a directory to sort
original_cifs = os.path.join(homeDirectory, "original_cifs")
createDirectory(original_cifs, "# WARN - No directory found for .cifs to sort. Place .cif files in or replace the newly created directory at", True)

cifs = {os.path.splitext(file)[0].replace(".cif", ""): file for file in sorted(os.listdir(original_cifs)) if file.endswith('.cif') and os.path.isfile(os.path.join(original_cifs, file))}

if len(cifs) == 0:#Make sure there are .cifs in the directory
    printToLog("# WARN - No .cif files found to sort. Place .cif files in ["+ original_cifs + "]")
    quit()
else:
    printToLog("# INFO - " + str(len(cifs)) + " .cif files found at ["+ original_cifs + "]")
    printToLog("# INFO - Following .cif files are available ["+str(list(cifs.values()))+"]")

#Make sure there is a directory for PSEUDOS
pseudosPath = os.path.join(homeDirectory, "PSEUDOS/")
createDirectory(pseudosPath, "# WARN - No directory found for PSEUDOS. Place .UPF files in or replace the newly created directory at", True)

psuedos = [file for file in os.listdir(pseudosPath) if file.endswith('.UPF') and os.path.isfile(os.path.join(pseudosPath, file))]#Get .UPFs from directory

if len(psuedos) == 0:#Make sure there are .UPFs in the directory
    printToLog("# WARN - No .UPF files found. Place .UPF files in ["+ pseudosPath + "]")
    quit()
else:
    printToLog("# INFO - " + str(len(psuedos)) + " .UPF files found at ["+ pseudosPath + "]")
    psuedo_elements = [file.split(".")[0] for file in psuedos]
    printToLog("# INFO - Following elements accounted for ["+str(list(set(psuedo_elements)))+"]")

dataAvailable = False

#Make sure there is a .csv to read from
structure_data = os.path.join(homeDirectory, "structure_data.csv")
if os.path.exists(structure_data):
    dataAvailable = True
    printToLog("# INFO - Loaded structure data from .csv at ["+ structure_data + "]")
else:
    printToLog("# INFO - No .csv file found to load compound data. Copy .csv from the CSD to the following path ["+ structure_data + "]")

printToLog("# INFO - Enter integer(s) with spaces between entries ('1 2 3') to choose processes to perform.")

OPT_R_VAL = "1"
OPT_DISORDER = "2"
OPT_SIZE = "3"
OPT_EMPTY = "4"
OPT_NO_H = "5"
OPT_INCOM_H = "6"
OPT_COCRYS = "7"

options = {
    OPT_R_VAL: f"Discard structures with r factor greater than [{rCap}]" if dataAvailable is not False else f"Not available. Provide structure_data.csv to filter by r value", 
    OPT_DISORDER: f"Discard structures flagged as disordered by their CSD author" if dataAvailable is not False else f"Not available. Provide structure_data.csv to filter by author flagged disorder",
    OPT_SIZE: f"Discard structures with volume greater than [{volumeCap}]",
    OPT_EMPTY: f"Discard structures without structural data",
    OPT_NO_H: f"Discard structures without hydrogen data",
    OPT_INCOM_H: f"Discard structures with incomplete hydrogen data",
    OPT_COCRYS: f"Discard structures with unreported cocrystals/solvent",
    "0": f"All in sequence",
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

removed = []
if not dataAvailable and choices.__contains__(OPT_R_VAL):
    choices.remove(OPT_R_VAL)
    removed.append(OPT_R_VAL)
if not dataAvailable and choices.__contains__(OPT_DISORDER):
    choices.remove(OPT_DISORDER)
    removed.append(OPT_DISORDER)
if len(removed) > 0:
    printToLog("# WARN - The following inputs ["+str(list(set(removed)))+"] are unavailable without a structure_data.csv")
    
printToLog("# INFO - The following processes have been selected ["+str(sorted(choices,key=int))+"]")

#Make sure there are directories to put the sorted files
cifs_path = os.path.join(homeDirectory, "cifs")
createDirectory(cifs_path, "# INFO - No directory found for storing .cifs, created at", False)

validated = os.path.join(cifs_path, "validated")
createDirectory(validated, "# INFO - No directory found for validated .cifs, created at", False)
validated_count = len(os.listdir(validated))

disordered = os.path.join(cifs_path, "disordered")
createDirectory(disordered, "# INFO - No directory found for disordered .cifs, created at", False)
disordered_count = len(os.listdir(disordered))

r_capped = os.path.join(cifs_path, "r_capped")
createDirectory(r_capped, "# INFO - No directory found for r value capped .cifs, created at", False)
r_capped_count = len(os.listdir(r_capped))

size_capped = os.path.join(cifs_path, "size_capped")
createDirectory(size_capped, "# INFO - No directory found for size capped .cifs, created at", False)
size_capped_count = len(os.listdir(size_capped))

unreported_cocrystal = os.path.join(cifs_path, "unreported_cocrystal")
createDirectory(unreported_cocrystal, "# INFO - No directory found for structure-formula mismatched .cifs, created at", False)
unreported_cocrystal_count = len(os.listdir(unreported_cocrystal))

empty = os.path.join(cifs_path, "empty")
createDirectory(empty, "# INFO - No directory found for empty .cifs, created at", False)
empty_count = len(os.listdir(empty))

no_hydrogen = os.path.join(cifs_path, "no_hydrogen")
createDirectory(no_hydrogen, "# INFO - No directory found for .cifs with no hydrogen data, created at", False)
no_hydrogen_count = len(os.listdir(no_hydrogen))

incomplete_hydrogen = os.path.join(cifs_path, "incomplete_hydrogen")
createDirectory(incomplete_hydrogen, "# INFO - No directory found for .cifs with incomplete hydrogen data, created at", False)
incomplete_hydrogen_count = len(os.listdir(incomplete_hydrogen))

unaccounted_elements = os.path.join(cifs_path, "unaccounted_elements")
createDirectory(unaccounted_elements, "# INFO - No directory found for .cifs with unaccounted element types, created at", False)
unaccounted_elements_count = len(os.listdir(unaccounted_elements))

if dataAvailable and choices.__contains__(OPT_R_VAL) or choices.__contains__(OPT_DISORDER):
    df = pd.read_csv(structure_data, encoding="utf-8-sig")
    df.set_index('[REFCODE]', inplace = True)


unaccountedElements = {}
for refcode, filename in cifs.items():
    unaccounted = []
    with open(os.path.join(original_cifs, filename), "r") as cif:
        printToLog("# INFO - Compound [" + refcode + "] Processing .cif file")
        formula_dict = {}
        counts = {}
        
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
                elif "_chemical_formula_sum" in line:
                    formula = re.match(r"_?(\w+)\s+(.*)", line).group(2).strip("\'")
                    printToLog("# INFO - Compound [" + refcode + "] Found to have formula ["+ formula +"]")

                    regex = re.compile('[^a-zA-Z ]')
                    atoms = list(set([atom.lower().capitalize() for atom in regex.sub('', formula).strip().split(" ")]))

                    printToLog("# INFO - Compound [" + refcode + "] Found to contain the following atom types: ["+ str(atoms) +"]")
                    for atom in atoms:
                        if not psuedo_elements.__contains__(atom):
                            unaccounted.append(atom)
                            if not atom in unaccountedElements:
                                unaccountedElements[atom] = 0
                            unaccountedElements[atom] += 1
                    
                    for element in formula.split(): #Save formula to a dict
                        formula_dict[re.findall('\\d+|\\D+', element)[0]] = int(re.findall('\\d+|\\D+', element)[1])
                    counts = formula_dict.copy()
                else:
                    if m := re.match(r"_?(\w+)\s+(.*)", line.strip()):
                        key = m.group(1)
                        value = m.group(2)
                        element_symbol = re.findall('\\d+|\\D+', value)[0].strip(" -")
                        if not "?" in key and counts.__contains__(element_symbol):
                            #Increment a copy of the formula dict for every corresponding atom position, counts down for convenience
                            counts[element_symbol] -= 1

    if len(unaccounted) > 0:
        printToLog("# INFO - Compound [" + refcode + "] discarded. Found to contain the following atoms not accounted for by .UPF files: ["+ str(unaccounted) +"]")
        shutil.copyfile(os.path.join(original_cifs, filename), os.path.join(unaccounted_elements, refcode + ".cif"))
        continue

    if choices.__contains__(OPT_R_VAL) or choices.__contains__(OPT_DISORDER):
        printToLog("# INFO - Compound ["+ refcode +"] Looking for structure_data.csv entry")
        if refcode in df.index:
            if choices.__contains__(OPT_R_VAL):
                rFactor = float(df.at[refcode, "[_refine_ls_R_factor]"])
                if rFactor > rCap:
                    printToLog("# INFO - Compound [" + refcode + "] Discarded. R factor [" + str(rFactor) + "] greater than cutoff point [" + str(rCap) + "]")
                    shutil.copyfile(os.path.join(original_cifs, filename), os.path.join(r_capped, refcode + ".cif"))
                    continue
            if choices.__contains__(OPT_DISORDER):
                disorder = str(df.at[refcode, "[_exptl_DISORDER]"])
                if not disorder == "nan":
                    printToLog("# INFO - Compound [" + refcode + "] Discarded. Disordered [" + str(disorder) + "]")
                    shutil.copyfile(os.path.join(original_cifs, filename), os.path.join(disordered, refcode + ".cif"))
                    continue
        else:
            printToLog("# WARN - Compound [" + refcode + "] Not present in structure_data.csv")
            continue  
    
    if choices.__contains__(OPT_SIZE):        
        if volume > volumeCap:
            printToLog("# INFO - Compound [" + refcode + "] Discarded. Cell Volume [" + str(volume) + "] greater than volume cap [" + str(volumeCap) + "]")
            shutil.copyfile(os.path.join(original_cifs, filename), os.path.join(size_capped, refcode + ".cif"))
            continue
    if not all(counts.get(key) % formula_dict.get(key) == 0 and not counts.get(key) == formula_dict.get(key) for key in formula_dict):
        if all(counts.get(key) == formula_dict.get(key) for key in formula_dict):
            if choices.__contains__(OPT_EMPTY):        
                printToLog("# INFO - Compound [" + refcode + "] Discarded. No structural data, all atoms unaccounted for: ["+ str(counts) +"]")
                shutil.copyfile(os.path.join(original_cifs, filename), os.path.join(empty, refcode + ".cif"))
                continue
        elif (counts.get("H") == formula_dict.get("H")): 
            if choices.__contains__(OPT_NO_H):        
                printToLog("# INFO - Compound [" + refcode + "] Discarded. No hydrogen data, unaccounted for atoms: ["+ str(counts) +"]")
                shutil.copyfile(os.path.join(original_cifs, filename), os.path.join(no_hydrogen, refcode + ".cif"))
                continue
        elif all(counts.get(key) == 0 or key == "H" for key in formula_dict):
            if choices.__contains__(OPT_INCOM_H):        
                printToLog("# INFO - Compound [" + refcode + "] Discarded. Incomplete hydrogen data, unaccounted for atoms: ["+ str(counts) +"]")
                shutil.copyfile(os.path.join(original_cifs, filename), os.path.join(incomplete_hydrogen, refcode + ".cif"))
                continue
        else:
            if choices.__contains__(OPT_COCRYS):        
                printToLog("# INFO - Compound [" + refcode + "] Discarded. Unreported co-crytsal/solvent, unaccounted for atoms: ["+ str(counts) +"]")
                shutil.copyfile(os.path.join(original_cifs, filename), os.path.join(unreported_cocrystal, refcode + ".cif"))
                continue
    
    printToLog("# INFO - Compound [" + refcode + "] Validated with constraints ["+str(sorted(choices,key=int))+"]")
    shutil.copyfile(os.path.join(original_cifs, filename), os.path.join(validated, refcode + ".cif"))


final = os.path.join(cifs_path, "original_cifs")

printToLog("# INFO - Moving sorted path ["+str(original_cifs)+"] to ["+str(final)+"]")
shutil.copytree(original_cifs, final, dirs_exist_ok=True)
removeDirectory(original_cifs, "# INFO - Cleaning sorted path at")

if choices.__contains__(OPT_R_VAL):
    printToLog(f"# INFO - [{'{: >3}'.format(len(os.listdir(r_capped)) - r_capped_count)}] compounds with r factor greater than [{rCap}]")
    
if choices.__contains__(OPT_DISORDER):
    printToLog(f"# INFO - [{'{: >3}'.format(len(os.listdir(disordered)) - disordered_count)}] compounds flagged as disordered by their CSD author")
    
if choices.__contains__(OPT_SIZE):
    printToLog(f"# INFO - [{'{: >3}'.format(len(os.listdir(size_capped)) - size_capped_count)}] compounds with volume greater than [{volumeCap}]")
    
if choices.__contains__(OPT_EMPTY):
    printToLog(f"# INFO - [{'{: >3}'.format(len(os.listdir(empty)) - empty_count)}] compounds without structural data")

if choices.__contains__(OPT_NO_H):
    printToLog(f"# INFO - [{'{: >3}'.format(len(os.listdir(no_hydrogen)) - no_hydrogen_count)}] compounds without hydrogen data")

if choices.__contains__(OPT_INCOM_H):
    printToLog(f"# INFO - [{'{: >3}'.format(len(os.listdir(incomplete_hydrogen)) - incomplete_hydrogen_count)}] compounds with incomplete hydrogen data")

if choices.__contains__(OPT_COCRYS):
    printToLog(f"# INFO - [{'{: >3}'.format(len(os.listdir(unreported_cocrystal)) - unreported_cocrystal_count)}] compounds with unreported cocrystals/solvent")

if len(unaccountedElements) > 0:
    printToLog(f"# INFO - [{'{: >3}'.format(len(os.listdir(unaccounted_elements)) - unaccounted_elements_count)}] compounds with elements unaccounted for by the current PSEUDOS")
    printToLog("# WARN - The following atoms are unaccounted for by the current PSUEDOS: ["+str(unaccountedElements)+"]")
else:
    printToLog("# INFO - All present atoms are accounted for by the current PSEUDOS")

surviving = len(os.listdir(validated)) - validated_count
printToLog("# INFO - Sort complete, [" + str(surviving) + "] of [" + str(len(cifs)) + "] compounds validated, with a % survival of [" +str(round((surviving / len(cifs)) * 100, 3))+ "%]")
