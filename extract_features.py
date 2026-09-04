import re 
import os
import sys
import time
import datetime
import pandas as pd
import numpy as np
import datetime
import time
import shutil
import math

#Jank thing to fix the path. Very annoying artefact of running python scripts by absolute path.
sys.path[0] = sys.path[0][:-6] + sys.path[0][-6:].replace("/utils", "")

from utils.generic_utils import printToLog as pl, createDirectory as cd, cellVolume, featureExtractor
from utils.params import *

from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor, rdFMCS, rdDistGeom
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeProperties, MMFFGetMoleculeForceField
from rdkit.Chem.rdMolTransforms import *

import octadist as oc

from utils.generic_utils import printToLog as pl, createDirectory as cd, getQueueLength, cellVolume, mol2Creator, featureExtractor

# Params
set_id = "MAIN" # set_id of the calculations to extract features from
pad = ["0","0","0","0"] # block to pad empty bond positions with (X atomic number, X bonds, W-P-X angle, P-X dist)

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

utils = os.path.join(homeDirectory,"utils")
atom_data = os.path.join(os.path.join(utils, "data"), "atom_data.csv")
if os.path.exists(atom_data):
    atom_data_df = pd.read_csv(atom_data, encoding="utf-8-sig")
    atom_data_df.set_index('Symbol', inplace = True)
else:
    printToLog("# WARN - No .csv file found to load atom data.")
    quit()

printToLog("# INFO - Enter integer(s) with spaces between entries ('1 2') to choose processes to perform.")
options = {
    "1": "Extract feature data from calculations",
    "2": "Extract feature data from experimental",
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
    
for choice in choices:    
    if not options.__contains__(choice):
        invalidInputs.append(choice)
if len(invalidInputs) > 0:
    printToLog("# WARN - The following inputs ["+str(list(set(invalidInputs)))+"] are not supported")
    quit()
printToLog("# INFO - The following processes have been selected ["+str(sorted(choices,key=int))+"]")

if choices.__contains__("1"):
    feature_data = os.path.join(homeDirectory, "_training_feature_data.csv")
    printToLog("# INFO - Creating new sheet ["+ feature_data + "]")

    feature_names = os.path.join(os.path.join(os.path.join(homeDirectory,"utils"), "data"), "feature_names.csv")
    if os.path.exists(feature_names):
        df = pd.read_csv(feature_names, encoding="utf-8-sig")        
        if os.path.isfile(feature_data):
            os.remove(feature_data)
        df.to_csv(feature_data, index=False)    
    else:
        printToLog("# WARN - No .csv file found to load feature names.")
        quit()
    
    output_files = os.path.join(os.path.join(homeDirectory, set_id), "output_files")
    createDirectory(output_files, "# WARN - No directory found for output files.", True)
    directories = sorted([directory for directory in os.listdir(output_files) if os.path.isdir(os.path.join(output_files, directory))])
    
    if len(directories) == 0:
        printToLog("# WARN - No directories found in ["+ output_files + "]")
        quit()
    else:
        printToLog(f"# INFO - Following directories are available [{directories}]")
        for refcode in directories:
            printToLog(f"# INFO - Processing compound [{refcode}]")
            
            refcodeDirectory = os.path.join(output_files, refcode)
            summary = os.path.join(refcodeDirectory, refcode+"_summary.txt")
            
            if os.path.exists(summary):
                with open(summary) as file:
                    summary_lines = file.readlines()
                    
                    join = ''.join(summary_lines)
                    if "BATCH_done = True" in join and "PWSCF_done = True" in join and "PWSCF_finalEnergy" in join and "GIPAW_done = True" in join:        
                        summary_atoms = join.split('#BEGIN_ATOMIC_POSITIONS\n')[1].split('\n#END_ATOMIC_POSITIONS')[0].split("\n")
                        for summary_number, summary_line in enumerate(summary_atoms, 0):
                            if summary_line.split()[1] == "W":
                                extractor = featureExtractor(log=log,directory=refcodeDirectory,refcode=refcode,summary_atoms=summary_atoms,site_id=summary_line.split()[0])
                                data = extractor.extract()
                                if not data == None: 
                                    df = pd.read_csv(feature_data)
                                                
                                    values = list(data.values())
                                    while len(values) < len(df.columns):
                                        values.extend(pad)
                                    df.loc[len(df)] = values

                                    df.to_csv(feature_data, index=False)    
                    else:
                        printToLog(f"# WARN - Compound [{refcode}] Not complete")
            else:
                printToLog(f"# WARN - Compound [{refcode}] No summary file found")       


    #Make sure there is a directory to sort
if choices.__contains__("2"):
    experimental = os.path.join(homeDirectory, "experimental")
    createDirectory(experimental, "# WARN - No directory found for experimental .cifs to process. Place .cif files in or replace the newly created directory at", True)
    
    feature_data = os.path.join(homeDirectory, "_experimental_feature_data.csv")                     
    feature_names = os.path.join(os.path.join(os.path.join(homeDirectory,"utils"), "data"), "feature_names.csv")

    if os.path.exists(feature_names):
        df = pd.read_csv(feature_names, encoding="utf-8-sig")        
        if os.path.isfile(feature_data):
            os.remove(feature_data)
        df.to_csv(feature_data, index=False)    
    else:
        printToLog("# WARN - No .csv file found to load feature names.")
        quit()
    
    
    cifs = {os.path.splitext(file)[0].replace(".cif", ""): file for file in sorted(os.listdir(experimental)) if file.endswith('.cif') and os.path.isfile(os.path.join(experimental, file))}
    
    if len(cifs) == 0:#Make sure there are .cifs in the directory
        printToLog("# WARN - No .cif files found to sort. Place .cif files in ["+ experimental + "]")
    else:
        printToLog("# INFO - " + str(len(cifs)) + " .cif files found at ["+ experimental + "]")
        printToLog("# INFO - Following .cif files are available ["+str(list(cifs.values()))+"]")
    
    for refcode, filename in cifs.items():
        refcodeDirectory = os.path.join(experimental, refcode)
        createDirectory(refcodeDirectory, f"# INFO - No directory found for compound [{refcode}], created at", False)
        shutil.copy(os.path.join(experimental, filename), os.path.join(refcodeDirectory, filename))
    
        cif = os.path.join(refcodeDirectory, filename)
    
        printToLog("# INFO - Compound [" + refcode + "] Sucessfully ran cif2cell")
        with open(cif) as file:    
            lines = file.readlines()
            
            cell_params = {}
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
                        cell_params[key] = str(round(value,4))
                    printToLog("# INFO - Compound ["+refcode+"] Cell params ["+str(cell_params)+"]")
        
            atom_positions = []
            for line in ''.join(lines).split('_atom_site_label\n')[1].split('\nloop_')[0].split("\n"):
                if not line.lstrip().startswith("_") and not line.lstrip().startswith("#") and not len(line.lstrip()) <= 9:
                    temp = re.sub('\s{2,}', ' ', line).lstrip().split()

                    printToLog(f"# INFO - Compound [{refcode}] found atom [{' '.join(temp[1:])}] from .cif")
                    atom_positions.append(' '.join(temp[1:]))
            
            creator = mol2Creator(log=log,directory=refcodeDirectory,refcode=refcode,cell_params=cell_params,atom_positions=atom_positions,df=atom_data_df)
            creator.create()
        
            for number, line in enumerate(atom_positions, 1):
                if line.split()[0].lower().capitalize() == "W":        
                    extractor = featureExtractor(log=log,directory=refcodeDirectory,refcode=refcode,summary_atoms=[""],site_id=f"#[{number}]")
                    data = extractor.extract()
                    if not data == None: 
                        df = pd.read_csv(feature_data)
                                    
                        values = list(data.values())
                        while len(values) < len(df.columns):
                            values.extend(pad)
                        df.loc[len(df)] = values
                        df.to_csv(feature_data, index=False)    



