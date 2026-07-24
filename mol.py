# mol.py - Jacob Duddridge

# Converts quantum-espresso output files into a supercell .cif file, preserving optimised positions
# converts supercell .cif to .mol2 using openbabel, automatically fixes broken bonds near metallic atoms
# extracts a single unit cell from .mol2, preserving connectivity

# All processes are reported to mol.log for debugging

#Imports
import os
import subprocess
import sys
import re 
import pandas as pd
import io
import csv
import datetime
import time
from utils.generic_utils import printToLog as pl, createDirectory as cd, writeCSV, getModules, isQueued
import shutil
import numpy as np
import math

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)

def append_atom(atom_number):
    atom_id = str('{: >6}'.format(atom_number))
    
    printToLog("# INFO - Compound ["+refcode+"] Atom appended ["+str(atoms[int(atom_number)].rstrip())+"]")
    new_atoms.append(atoms[int(atom_number)])
        
    for bond_number, bond_line in enumerate(bonds, 1):
        first_id = str(bond_line[6:12])                                        
        second_id = str(bond_line[12:18])

        if atom_id == first_id:
            if not new_bonds.__contains__(bond_line):
                new_bonds.append(bond_line)

            second_id = int(second_id.lstrip())
            if not new_atoms.__contains__(atoms[second_id]):
                append_atom(second_id)
            else:
                printToLog("# INFO - Compound ["+refcode+"] Atom already selected ["+str(atoms[second_id].rstrip())+"]")

        if atom_id == second_id:
            if not new_bonds.__contains__(bond_line):
                new_bonds.append(bond_line)

            first_id = int(first_id.lstrip())
            if not new_atoms.__contains__(atoms[first_id]):
                append_atom(first_id)
            else:
                printToLog("# INFO - Compound ["+refcode+"] Atom already selected ["+str(atoms[first_id].rstrip())+"]")

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

utils = os.path.join(homeDirectory,"utils")
atom_data = os.path.join(utils, "atom_data.csv")
if os.path.exists(atom_data):
    df = pd.read_csv(atom_data, encoding="utf-8-sig")
    df.set_index('Symbol', inplace = True)
else:
    printToLog("# WARN - No .csv file found to load atom data.")
    quit()
    
#Make sure there is a directory to process
inputPath = os.path.join(homeDirectory, "Output_Files")
createDirectory(inputPath, "# WARN - No directory found for output files.", True)
directories = [directory for directory in os.listdir(inputPath) if os.path.isdir(os.path.join(inputPath, directory)) and not directory.startswith(".") and os.path.isfile(os.path.join(os.path.join(inputPath, directory), directory+".out")) and os.path.isfile(os.path.join(os.path.join(inputPath, directory), "gipaw."+directory+".out"))]

numberOfDirectories = len(directories) # determine number of directories
if numberOfDirectories == 0:
    printToLog("# WARN - No directories found in ["+ inputPath + "]")
    quit()
else:
    printToLog("# INFO - [" + str(numberOfDirectories) + "] directories found at ["+ inputPath + "]")
    directories = sorted(directories)
    for refcode in directories:
        printToLog("# INFO - Processing compound with refcode ["+ refcode +"]")
        if not isQueued(log, refcode):
            refcodeDirectory = os.path.join(inputPath, refcode)
        
            out = os.path.join(refcodeDirectory, refcode+".out")    
            cif = os.path.join(refcodeDirectory, refcode+".cif")
            opt_cif = os.path.join(refcodeDirectory, refcode+"_opt.cif")
            super_mol2 = os.path.join(refcodeDirectory, refcode+"_super.mol2")
            
            mol2 = os.path.join(refcodeDirectory, refcode+".mol2")
            mol2_copy = os.path.join(os.path.join(homeDirectory, "mol"), refcode+".mol2")                    
            
            if os.path.exists(opt_cif):
                printToLog("# INFO - Compound ["+ refcode +"] Cleaning existing optimised _super.cif file ["+ opt_cif + "]")
                os.remove(opt_cif)
            if os.path.exists(super_mol2):
                printToLog("# INFO - Compound ["+ refcode +"] Cleaning existing _super.mol2 file ["+ super_mol2 + "]")
                os.remove(super_mol2)
            if os.path.exists(mol2):
                printToLog("# INFO - Compound ["+ refcode +"] Cleaning existing .mol2 file ["+ mol2 + "]")
                os.remove(mol2)
            if os.path.exists(mol2_copy):
                printToLog("# INFO - Compound ["+ refcode +"] Cleaning existing .mol2 file copy ["+ mol2_copy + "]")
                os.remove(mol2_copy)

            single_cell_count = 0
            if os.path.isfile(cif):
                if os.path.isfile(out):
                    with open(opt_cif, "w") as opt:
                        printToLog("# INFO - Compound ["+refcode+"] Populating optimised .cif file")
                        print("data_"+refcode+"_OPT", file=opt)
    
                        with open(cif) as file:
                            lines = file.readlines()
                            for number, line in enumerate(lines, 0):   
                                if "_chemical_formula_sum" in line:
                                    print(line.strip(), file=opt)
                                elif "_cell_length_a" in line:
                                    print(lines[number].strip(), file=opt)
                                    print(lines[number+1].strip(), file=opt)
                                    print(lines[number+2].strip(), file=opt)
                                    print(lines[number+3].strip(), file=opt)
                                    print(lines[number+4].strip(), file=opt)
                                    print(lines[number+5].strip(), file=opt)
    
                        print("_space_group_name_H-M_alt 'P 1'", file=opt)
                        print("_space_group_IT_number 1", file=opt)
                        print("loop_", file=opt)
                        print("_space_group_symop_operation_xyz", file=opt)
                        print("'x, y, z'", file=opt)
    
                        print("loop_", file=opt)
                        print("_atom_site_label", file=opt)
                        print("_atom_site_type_symbol", file=opt)
                        print("_atom_site_fract_x", file=opt)
                        print("_atom_site_fract_y", file=opt)
                        print("_atom_site_fract_z", file=opt)

                        start = 0
                        end = 0
                        
                        printToLog("# INFO - Compound ["+refcode+"] Processing .out file")
                        with open(out) as file:                
                            lines = file.readlines()
                            for number, line in enumerate(lines, 1):      
                                if "Begin final coordinates" in line:
                                    start = number + 2                    
                                elif "End final coordinates" in line:
                                    end = number - 1

                            single_cell_count = len(lines[start:end])
                            if single_cell_count == 0:
                                printToLog("# WARN - Compound ["+refcode+"] Does not have final coordinates")
                                continue
                            
                            counts = {}

                            arrays = []
                            arrays.append([0, 0, 0])

                            arrays.append([1, 1, 1])
                            arrays.append([1, 1, 0])
                            arrays.append([1, 1, -1])
                            arrays.append([1, 0, 1])
                            arrays.append([1, 0, 0])
                            arrays.append([1, 0, -1])
                            arrays.append([1, -1, 1])
                            arrays.append([1, -1, 0])
                            arrays.append([1, -1, -1])
                            arrays.append([0, 1, 1])
                            arrays.append([0, 1, 0])
                            arrays.append([0, 1, -1])
                            arrays.append([0, 0, 1])
                           #arrays.append([0, 0, 0])
                            arrays.append([0, 0, -1])
                            arrays.append([0, -1, 1])
                            arrays.append([0, -1, 0])
                            arrays.append([0, -1, -1])
                            arrays.append([-1, 1, 1])
                            arrays.append([-1, 1, 0])
                            arrays.append([-1, 1, -1])
                            arrays.append([-1, 0, 1])
                            arrays.append([-1, 0, 0])
                            arrays.append([-1, 0, -1])
                            arrays.append([-1, -1, 1])
                            arrays.append([-1, -1, 0])
                            arrays.append([-1, -1, -1])

                            printToLog("# INFO - Compound ["+refcode+"] Creating _super.cif")
                            print("#ATOMS_START", file=opt)
                            for array in arrays:
                                print("#START"+str(array), file=opt)
                                for number, line in enumerate(lines[start:end], 0):
                                    curr = re.sub('\s{2,}', ' ', line).split()
                                    element = curr[0].lower().capitalize()
    
                                    if not counts.__contains__(element):
                                        counts[element] = 0
                                    counts[element] += 1
                                    new = element + str(counts[element]) + " " + element + " " + str(float(curr[1])+array[0]) + " " + str(float(curr[2])+array[1]) + " " + str(float(curr[3])+array[2])
                                    print(new, file=opt)
                                print("#END"+str(array), file=opt)                        
                            print("#ATOMS_END", file=opt)

                    printToLog("# INFO - Compound ["+refcode+"] Attempting to create _super.mol2")
                    try:
                        subprocess.call(f"module load {getModules()}; cd {refcodeDirectory}; obabel -i cif {refcode}_opt.cif -o mol2 -O {refcode}_super.mol2",shell=True)
                        printToLog("# INFO - Compound ["+refcode+"] Successfully created _super.mol2")

                        if os.path.isfile(super_mol2):
                            with open(super_mol2, "r") as opt:
                                lines = opt.readlines()
                            atoms_positions = ''.join(lines).split('@<TRIPOS>ATOM')[1].split('@<TRIPOS>BOND')[0]

                            count = 0
                            with open(super_mol2, "w") as file:
                                for line in lines:
                                    line = line.rstrip("\n")
                                    if line in atoms_positions and not line == "":
                                        count += 1
                                        if count > single_cell_count:
                                            count = 1
                                        line += " #" + str(count)
                                    print(line, file=file)


                            with open(super_mol2) as file:
                                lines = file.readlines()
                            with open(super_mol2, "a") as file:
                                atom_line = 0
                                bond_line = 0
                                for number, line in enumerate(lines, 1):                  
                                    if "@<TRIPOS>ATOM" in line:
                                        atom_line = number
                                    if "@<TRIPOS>BOND" in line:
                                        bond_line = number

                                bonds = lines[bond_line:]
                                atoms = lines[atom_line:bond_line-1]
                
                                bond_count = len(bonds)
                                printToLog("# INFO - Compound ["+refcode+"] Attempting to find metallic bonds")
                                for metal_number, metal_line in enumerate(atoms, 1):
                                    metal_curr = re.sub('\s{2,}', ' ', metal_line).split()
                                    metal = metal_curr[1].lower().capitalize()
                
                                    if str(df.at[metal, "Metal"]) == "True":
                                        metal_array = np.array([float(metal_curr[2]), float(metal_curr[3]), float(metal_curr[4])])
                                        for other_number, other_line in enumerate(atoms, 1):
                                            other_curr = re.sub('\s{2,}', ' ', other_line).split()
                                            other = other_curr[1].lower().capitalize()
                                            
                                            if not other == "H" and not metal_number == other_number: 
                                                other_array = np.array([float(other_curr[2]), float(other_curr[3]), float(other_curr[4])])
                                                cutoff = float(df.at[metal, "RCov"]) + float(df.at[other, "RCov"]) + 0.45
                                                distance = np.sqrt(np.sum((metal_array-other_array)**2, axis=0))
                                                if distance < cutoff:
                                                    metal_id = str('{: >6}'.format(metal_number))
                                                    other_id = str('{: >6}'.format(other_number))
                
                                                    bond_exists = False
                                                    for bond_number, bond_line in enumerate(bonds, 1):
                                                        existing_ids = str(bond_line[6:18])
                                                        if existing_ids == metal_id+other_id or existing_ids == other_id+metal_id:
                                                            bond_exists = True
                                                    if not bond_exists:
                                                        bond_count += 1
                
                                                        printToLog("# INFO - Compound ["+refcode+"] Adding bond ["+ str('{: >6}'.format(metal_number))+" "+ str(metal)+" - "+ str('{: >6}'.format(other_number))+" "+ str(other)+ " "+str(round(cutoff,3))+" "+ str(round(distance,3))+"]")
                                                        print(str('{: >6}'.format(bond_count)) + str('{: >6}'.format(metal_number)) + str('{: >6}'.format(other_number)) + "    1", file=file)
                                                    else:
                                                        printToLog("# INFO - Compound ["+refcode+"] Bond already present ["+ str('{: >6}'.format(metal_number))+" "+ str(metal)+" - "+ str('{: >6}'.format(other_number))+" "+ str(other)+ " "+str(round(cutoff,3))+" "+ str(round(distance,3))+"]")

                            printToLog("# INFO - Compound ["+refcode+"] Added ["+str(bond_count - len(bonds))+"] missing metallic bonds")
                            printToLog("# INFO - Compound ["+refcode+"] Exploring bonding networks")
                         
                            atom_start_line = 0
                            bond_start_line = 0
                            with open(mol2, "w") as mol2_file:
                                new_atoms = []
                                new_bonds = []
                                
                                with open(super_mol2) as file:
                                    lines = file.readlines()

                                    for number, line in enumerate(lines, 1):                  
                                        if "@<TRIPOS>ATOM" in line:
                                            atom_start_line = number
                                        if "@<TRIPOS>BOND" in line:
                                            bond_start_line = number
                                    
                                    bonds = lines[bond_start_line:]
                                    atoms = lines[atom_start_line-1:bond_start_line-1]
                                    for number, line in enumerate(lines[:atom_start_line], 0):
                                        print(line.rstrip(), file=mol2_file)

                                    for atom_number, atom_line in enumerate(lines[atom_start_line:int(atom_start_line+single_cell_count)], 1):
                                        if not new_atoms.__contains__(atom_line):
                                            append_atom(atom_number)
                                        else:
                                            printToLog("# INFO - Compound ["+refcode+"] Atom already selected ["+str(atoms[atom_number].rstrip())+"]")
                                
                                printToLog("# INFO - Compound ["+refcode+"] Populating .mol2")
                                id_map = {}
                                for atom_number, new_atom in enumerate(new_atoms, 1):
                                    new_id = str('{: >6}'.format(atom_number))
                                    old_id = str(new_atom[1:7])

                                    id_map[old_id] = new_id
                                    temp = " " + str(new_id) + new_atom.rstrip()[7:]
                                    print(temp, file=mol2_file)
                                    
                                print("@<TRIPOS>BOND", file=mol2_file)
                                for bond_number, new_bond in enumerate(new_bonds, 1):
                                    first_id = str(new_bond[6:12])                                        
                                    second_id = str(new_bond[12:18])       
                                    
                                    new_id = str('{: >6}'.format(bond_number))
                                    temp = new_id + id_map[first_id] + id_map[second_id] + new_bond[18:].rstrip()
                                    print(temp, file=mol2_file)

                                printToLog("# INFO - Compound ["+ refcode +"] Contains ["+str(len(new_atoms))+"] atoms and ["+str(len(new_bonds))+"] bonds")
                                if (len(new_bonds) / len(new_atoms)) < 0.8:
                                    printToLog("# WARN - Compound ["+ refcode +"] Has unusually low bond density")
                                
                            with open(mol2) as file:
                                lines = file.readlines()
                            with open(mol2, "w") as file:
                                for number, line in enumerate(lines, 1):
                                    if number == 3: #MAYBE BETTER ANSWER FOR THIS?
                                        temp = line.strip().lstrip().split()                                        
                                        print(" "+str(len(new_atoms))+" "+str(len(new_bonds)) +" "+str(temp[2])+" "+str(temp[3])+" "+str(temp[4]), file=file)
                                    else:
                                        print(line.rstrip(), file=file)
                                        
                            shutil.copyfile(mol2, mol2_copy)
                            printToLog("# INFO - Compound ["+ refcode +"] Copied .mol2 file ["+ mol2_copy + "]")
                        else:
                            printToLog("# WARN - Compound ["+refcode+"] .mol2 output not found")
                    except subprocess.CalledProcessError as e:
                        printToLog("# WARN - Compound ["+refcode+"] Error creating _super.mol2")
                        printToLog(str(e))   
                else:
                    printToLog("# WARN - Compound ["+refcode+"] No .out file found")
            else:
                printToLog("# WARN - Compound ["+refcode+"] No .cif file found")

            
