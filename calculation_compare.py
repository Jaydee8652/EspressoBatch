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

tolerance = 0.01
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
qe_params = os.path.join(homeDirectory, "qe_params.csv")

if not os.path.isfile(qe_params):
    printToLog(log,"# WARN - No qe_params.csv found.")
    quit()

#Make sure there is a directory to process
cifs_path = os.path.join(homeDirectory, "cifs")
createDirectory(cifs_path, "# WARN - No directory found for .cifs to process.", False)

validated = os.path.join(cifs_path, "validated")
createDirectory(validated, "# WARN - No directory found for .cifs to process. Place .cif files in or replace the newly created directory at ["+validated+"]", True)

cifs = [os.path.splitext(file)[0].replace(".cif", "") for file in sorted(os.listdir(validated)) if file.endswith('.cif') and os.path.isfile(os.path.join(validated, file))]

qe = pd.read_csv(qe_params, encoding="utf-8-sig")
qe.set_index('set_id', inplace = True)

for refcode in cifs:
    printToLog(f"# INFO - Processing [{refcode}]")

    global_atoms = {}
    global_atom_ids = []
    for set_id, row in qe.iterrows():
        printToLog(f"# INFO - Processing set_id [{set_id}]")
    
        id_directory = os.path.join(homeDirectory, set_id)
        summary_file = os.path.join(os.path.join(id_directory, "summary_files"), f"{refcode}_summary.txt")
        
        if os.path.isfile(summary_file):
            with open(summary_file) as file:
                lines = file.readlines()
                atom_lines = ''.join(lines).split('#BEGIN_ATOMIC_POSITIONS\n')[1].split('\n#END_ATOMIC_POSITIONS')[0].split("\n")
                ids = list(map(lambda atom: ' '.join(atom.split()[0:2]), atom_lines))
                
                global_atoms[set_id] = dict(zip(ids, atom_lines))
                global_atom_ids.append(ids)                             
        
    if len(global_atom_ids) > 1 and global_atom_ids[:-1] == global_atom_ids[1:]:
        for atom_id in global_atom_ids[0]:
            for i in range(len(global_atoms.items())):
                for j in range(i + 1, len(global_atoms.items())):
                    main = list(global_atoms.items())[i][1][atom_id]
                    comp = list(global_atoms.items())[j][1][atom_id]
        
                    main_set_id = list(global_atoms.items())[i][0]
                    comp_set_id = list(global_atoms.items())[j][0]

                    main_sigma = float(main.split("(")[1].split(")")[0])
                    comp_sigma = float(comp.split("(")[1].split(")")[0])
                    sigma_diff = abs(main_sigma - comp_sigma)
                    if sigma_diff > tolerance:
                        printToLog(f"# WARN - Compound [{refcode}] has equivalent atoms [{atom_id}] with sigma values [{main_set_id} ({main_sigma}) {comp_set_id} ({comp_sigma})] outside tolerance [{sigma_diff}]")
                    

                
                #data[f"W sigma"] = atom.split("(")[1].split(")")[0]
                #data[f"W sigma_33"] = atom.split("[XYZ ")[1].split("]")[0].split()
                #data[f"W sigma_11"] = atom.split("[sigma_11 ")[1].split("]")[0]
                #data[f"W sigma_22"] = atom.split("[sigma_22 ")[1].split("]")[0]
                #data[f"W sigma_33"] = atom.split("[sigma_33 ")[1].split("]")[0]