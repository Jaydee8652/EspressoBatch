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

cifs_path = os.path.join(homeDirectory, "cifs")
createDirectory(cifs_path, "# WARN - No directory found for .cifs to process.", False)

validated = os.path.join(cifs_path, "validated")
createDirectory(validated, "# WARN - No directory found for .cifs to process. Place .cif files in or replace the newly created directory at ["+validated+"]", True)


cifs = {os.path.splitext(file)[0].replace(".cif", ""): file for file in sorted(os.listdir(validated)) if file.endswith('.cif') and os.path.isfile(os.path.join(validated, file))}

if len(cifs) == 0:#Make sure there are .cifs in the directory
    printToLog("# WARN - No .cif files found to sort. Place .cif files in ["+ validated + "]")
    quit()
else:
    printToLog("# INFO - " + str(len(cifs)) + " .cif files found at ["+ validated + "]")
    printToLog("# INFO - Following .cif files are available ["+str(list(cifs.values()))+"]")

#Make sure there is a .csv to read from
structure_data = os.path.join(homeDirectory, "structure_data.csv")
if os.path.exists(structure_data):
    printToLog("# INFO - Loaded structure data from .csv at ["+ structure_data + "]")
else:
    printToLog("# INFO - No .csv file found to load compound data. Copy .csv from the CSD to the following path ["+ structure_data + "]")
    quit()

df = pd.read_csv(structure_data, encoding="utf-8-sig")
df.set_index('[REFCODE]', inplace = True)

global_dict = {}
for refcode, filename in cifs.items():
    printToLog("# INFO - Compound ["+ refcode +"] Looking for structure_data.csv entry")
    if refcode in df.index:
        authors = str(df.at[refcode, "[_publ_authors]"]).split(",")
        journal = str(df.at[refcode, "[_journal_name_full]"]).strip()
        for author in authors:
            author = author.strip()
            if not global_dict.__contains__(author):
                global_dict[author] = {}
                global_dict[author]["Total"] = 0
            
            if not global_dict[author].__contains__(journal):
               global_dict[author][journal] = 0
            global_dict[author][journal] += 1 
            global_dict[author]["Total"] += 1 

with open(os.path.join(homeDirectory, "author_summary.txt"), "a") as file:
    for key in reversed(sorted(global_dict, key=lambda x: (global_dict[x]['Total']))):
        print(f"{key} {global_dict[key]}")
        print(f"{key} {global_dict[key]}", file=file)

