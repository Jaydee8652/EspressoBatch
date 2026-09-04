#Generic utility functions used by all scripts
import csv
import datetime
import time
import os
import shutil
import subprocess
import numpy as np
import re
import sys

#Jank thing to fix the path. Very annoying artefact of running python scripts by absolute path.
sys.path[0] = sys.path[0][:-6] + sys.path[0][-6:].replace("/utils", "")

from utils.params import * 

logs = os.path.join(os.getcwd(), "logs")
if not os.path.exists(logs):
    os.makedirs(logs)

#Prints and logs in one, convention I personally like
def printToLog(log, info):
    time = ""
    if not str(info).startswith(" ---"):
        time = str(datetime.datetime.now().strftime("[%H:%M:%S]"))
    
    with open(os.path.join(logs, log), "a") as log:
        print(f"{time} {info}",file=log)
        print(f"{time} {info}")
        
#Create directory if it doesn't exist. Optionally crash deliberately if doesn't exist
def createDirectory(log, path, text, exit):
    if not os.path.exists(path):
        printToLog(log, text + " ["+ path + "]")
        os.makedirs(path)
        if exit:
            quit()

#Remove directory if it exists
def removeDirectory(log, path, text):
    if os.path.exists(path):
        printToLog(log, text + " ["+ path + "]")
        shutil.rmtree(path)

#Write an entry to a local csv
def writeCSV(df, refcode, location, value):
    if not value == "":
        df.loc[refcode, location] = value

# Calculate cell volume from cell params
def cellVolume(cell_a, cell_b, cell_c, cell_α, cell_β, cell_γ):
    return cell_a * cell_b * cell_c * np.sqrt((1 - (np.cos(cell_α) ** 2) - (np.cos(cell_β) ** 2) - (np.cos(cell_γ)) ** 2) + (2 * np.cos(cell_α) * np.cos(cell_β) * np.cos(cell_γ)))

# Convert alat + vectors to list of cell params
# returns in order [cell_a, cell_b, cell_c, cell_α, cell_β, cell_γ]
def parseAlat(Alat, vector1, vector2, vector3):    
    a = np.multiply(vector1, Alat)
    b = np.multiply(vector2, Alat)
    c = np.multiply(vector3, Alat)

    return [np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c), vectorAngle(b,c), vectorAngle(a,c), vectorAngle(a,b)]

# Calculate angle between two 3D vectors
def vectorAngle(A, B):
    dot_product = np.dot(A, B)
    magnitude_A = np.linalg.norm(A)
    magnitude_B = np.linalg.norm(B)
    
    return np.degrees(np.arccos(dot_product / (magnitude_A * magnitude_B)))

#Get and decode the current slurm queue. Can be read like a file
def getQueue(log):
    printToLog(log,"# INFO - Attempting to retrieve current slurm queue.")
    try:
        # Uses custom output formatting so that job names (in theory) don't get truncated. Refcodes longer than 41 characters will break
        # Increase 45j if necessary
        out = subprocess.check_output(['squeue --format="%.10i %.10a %.45j %.2t %.10L %.10M %.6C %.6D %.6m %R" --me'],shell=True)
        out = out.decode("utf-8")
        return out
    except subprocess.CalledProcessError as e:
        printToLog(log,"# INFO - Error retreiving slurm queue.")
        printToLog(log,str(e))

# Only considers jobs ending in _SUB so that other jobs such as jupyter sessions don't interfere
def getQueued(log):
    lines = getQueue(log).splitlines()

    active_jobs = []
    for line in lines:
        printToLog(log, line)

        job_id = re.sub('\s{2,}', ' ', line).strip().split(" ")[2]
        if "_SUB" in job_id:
            active_jobs.append(job_id)
        
    return active_jobs

# Get the length of the current slurm queue
def getQueueLength(log):
    length = len(getQueued(log))
    printToLog(log,"# INFO - Slurm queue contains ["+str(length)+"] batched calculations.")
    return length

#Check if a specific refcode is in the queue
def isQueued(log, refcode):
    if getQueued(log).__contains__(refcode+"_SUB"):
        printToLog(log,"# INFO - Compound ["+refcode+"] is currently queued.")
        return True
    printToLog(log,"# INFO - Compound ["+refcode+"] is not currently queued.")
    return False

from rdkit import Chem
from rdkit.Chem.rdMolTransforms import *
import octadist as oc

class mol2Creator:
    #Functions
    def printToLog(self,info):#Prints and logs in one, convention I personally like
        printToLog(self.log, info)    

    def populate_mol(self, path, lines, seed=1):           
        with open(path, "w") as file:                        
            for number, line in enumerate(''.join(lines).split("\n@<TRIPOS>ATOM")[0].split("\n"), 0):
                print(line.rstrip(), file=file)
    
            if seed == 1:
                self.append_atom(atom_number=seed)
            else:   
                for atom_number, atom_line in enumerate(seed, 1):
                    if not self.new_atoms.__contains__(atom_line):
                        self.append_atom(atom_number=atom_number)
                    else:
                        self.printToLog("# INFO - Compound ["+self.refcode+"] Atom already selected ["+str(self.atoms[atom_number].rstrip())+"]")
    
            self.printToLog(f"# INFO - Compound [{self.refcode}] Populating single cell .mol2" if not seed == 1 else f"# INFO - Compound [{self.refcode}] Populating single molecule .mol2")
            id_map = {}
            
            print("@<TRIPOS>ATOM", file=file)
            for atom_number, new_atom in enumerate(self.new_atoms, 1):
                new_id = str('{: >6}'.format(atom_number))
                old_id = str(new_atom[1:7])
    
                id_map[old_id] = new_id
                temp = " " + str(new_id) + new_atom.rstrip()[7:]
                print(temp, file=file)
                
            print("@<TRIPOS>BOND", file=file)
            for bond_number, new_bond in enumerate(self.new_bonds, 1):
                first_id = str(new_bond[6:12])                                        
                second_id = str(new_bond[12:18])       
                
                new_id = str('{: >6}'.format(bond_number))
                temp = new_id + id_map[first_id] + id_map[second_id] + new_bond[18:].rstrip()
                print(temp, file=file)
    
            self.printToLog("# INFO - Compound ["+ self.refcode +"] Contains ["+str(len(self.new_atoms))+"] atoms and ["+str(len(self.new_bonds))+"] bonds")
            if (len(self.new_bonds) / len(self.new_atoms)) < 0.8:
                self.printToLog("# WARN - Compound ["+ self.refcode +"] Has unusually low bond density")
                
        with open(path) as file:
            lines = file.readlines()
        with open(path, "w") as file:
            for number, line in enumerate(lines, 1):
                if number == 3: #MAYBE BETTER ANSWER FOR THIS?
                    temp = line.strip().lstrip().split()                                        
                    print(" "+str(len(self.new_atoms))+" "+str(len(self.new_bonds)) +" "+str(temp[2])+" "+str(temp[3])+" "+str(temp[4]), file=file)
                else:
                    print(line.rstrip(), file=file)
    
    def append_atom(self, atom_number):
        atom_id = str('{: >6}'.format(atom_number))
        
        self.printToLog("# INFO - Compound ["+self.refcode+"] Atom appended ["+str(self.atoms[int(atom_number)-1].rstrip())+"]")
        self.new_atoms.append(self.atoms[int(atom_number)-1])
            
        for bond_number, bond_line in enumerate(self.bonds, 1):
            first_id = str(bond_line[6:12])                                        
            second_id = str(bond_line[12:18])
    
            if atom_id == first_id:
                if not self.new_bonds.__contains__(bond_line):
                    self.new_bonds.append(bond_line)
    
                second_id = int(second_id.lstrip())
                if not self.new_atoms.__contains__(self.atoms[second_id-1]):
                    self.append_atom(atom_number=second_id)
                else:
                    self.printToLog("# INFO - Compound ["+self.refcode+"] Atom already selected ["+str(self.atoms[second_id-1].rstrip())+"]")
    
            if atom_id == second_id:
                if not self.new_bonds.__contains__(bond_line):
                    self.new_bonds.append(bond_line)
    
                first_id = int(first_id.lstrip())
                if not self.new_atoms.__contains__(self.atoms[first_id-1]):
                    self.append_atom(atom_number=first_id)
                else:
                    self.printToLog("# INFO - Compound ["+self.refcode+"] Atom already selected ["+str(self.atoms[first_id-1].rstrip())+"]")
     
    def __init__(self, log, directory, refcode, cell_params, atom_positions, df):
        self.log = log
        self.directory = directory
        self.refcode = refcode
        self.cell_params = cell_params
        self.atom_positions = atom_positions
        self.df = df
   
    def create(self):
        opt_cif = os.path.join(self.directory, self.refcode+"_opt.cif")
        super_mol2 = os.path.join(self.directory, self.refcode+"_super.mol2")
        cell_mol2 = os.path.join(self.directory, self.refcode+"_cell.mol2")
        molecule_mol2 = os.path.join(self.directory, self.refcode+"_molecule.mol2")
    
        with open(opt_cif, "w") as opt:        
            self.printToLog("# INFO - Compound ["+self.refcode+"] Populating optimised .cif file")
            
            print("data_"+self.refcode+"_OPT", file=opt)
            print("_cell_length_a " + self.cell_params["a"],file=opt)
            print("_cell_length_b " + self.cell_params["b"],file=opt)
            print("_cell_length_c " + self.cell_params["c"],file=opt)
            print("_cell_angle_alpha " + self.cell_params["α"],file=opt)
            print("_cell_angle_beta " + self.cell_params["β"],file=opt)
            print("_cell_angle_gamma " + self.cell_params["γ"],file=opt)
            
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
    
            counts = {}
            print("#ATOMS_START", file=opt)
            for array in arrays:
                print("#START"+str(array), file=opt)
                for number, line in enumerate(self.atom_positions, 0):
                    curr = re.sub('\s{2,}', ' ', line).split()
                    element = curr[0].lower().capitalize()
    
                    if not counts.__contains__(element):
                        counts[element] = 0
                    counts[element] += 1
    
                    new = element + str(counts[element]) + " " + element + " " + str(float(curr[1].split("(")[0])+array[0]) + " " + str(float(curr[2].split("(")[0])+array[1]) + " " + str(float(curr[3].split("(")[0])+array[2])
                    print(new, file=opt)
                print("#END"+str(array), file=opt)                        
            print("#ATOMS_END", file=opt)
            self.printToLog("# INFO - Compound ["+self.refcode+"] Created _super.cif")
        
        try: 
            subprocess.call(f"module load {param_modules}; cd {self.directory}; obabel -i cif {self.refcode}_opt.cif -o mol2 -O {self.refcode}_super.mol2",shell=True)

            #subprocess.call(f"module load {param_modules}; cd {refcodeDirectory}; obabel -i cif {refcode}_opt.cif -o mol2 -O {refcode}_super.mol2",shell=True)

            self.printToLog("# INFO - Compound ["+self.refcode+"] Created _super.mol2")
        
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
                            if count > len(self.atom_positions):
                                count = 1
                            line += f" #[{count}]"
                        print(line, file=file)
        
        
                with open(super_mol2) as file:
                    lines = file.readlines()
                with open(super_mol2, "a") as file:
                    atoms = ''.join(lines).split('@<TRIPOS>ATOM\n')[1].split('\n@<TRIPOS>BOND')[0].split("\n")
                    bonds = ''.join(lines).split('@<TRIPOS>BOND\n')[1].split("\n")
        
                    bond_count = len(bonds)
                    self.printToLog("# INFO - Compound ["+self.refcode+"] Attempting to find metallic bonds")
                    for metal_number, metal_line in enumerate(atoms, 1):
                        metal_curr = re.sub('\s{2,}', ' ', metal_line).split()
                        metal = metal_curr[1].lower().capitalize()
        
                        if str(self.df.at[metal, "Metal"]) == "True":
                            metal_array = np.array([float(metal_curr[2]), float(metal_curr[3]), float(metal_curr[4])])
                            for other_number, other_line in enumerate(atoms, 1):
                                other_curr = re.sub('\s{2,}', ' ', other_line).split()
                                other = other_curr[1].lower().capitalize()
                                
                                if not other == "H" and not metal_number == other_number: 
                                    other_array = np.array([float(other_curr[2]), float(other_curr[3]), float(other_curr[4])])
                                    cutoff = float(self.df.at[metal, "RCov"]) + float(self.df.at[other, "RCov"]) + 0.45
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
        
                                            self.printToLog("# INFO - Compound ["+self.refcode+"] Adding bond ["+ str('{: >6}'.format(metal_number))+" "+ str(metal)+" - "+ str('{: >6}'.format(other_number))+" "+ str(other)+ " "+str(round(cutoff,3))+" "+ str(round(distance,3))+"]")
                                            
                                            bond = str('{: >6}'.format(bond_count)) + str('{: >6}'.format(metal_number)) + str('{: >6}'.format(other_number)) + "    1"
                                            print(bond, file=file)
                                            bonds.append(bond)
                                        else:
                                            self.printToLog("# INFO - Compound ["+self.refcode+"] Bond already present ["+ str('{: >6}'.format(metal_number))+" "+ str(metal)+" - "+ str('{: >6}'.format(other_number))+" "+ str(other)+ " "+str(round(cutoff,3))+" "+ str(round(distance,3))+"]")
        
                self.printToLog("# INFO - Compound ["+self.refcode+"] Added ["+str(bond_count - len(bonds))+"] missing metallic bonds")
                self.printToLog("# INFO - Compound ["+self.refcode+"] Exploring bonding networks")
        
                with open(super_mol2) as file:
                    lines = file.readlines()
                    self.atoms = ''.join(lines).split('@<TRIPOS>ATOM\n')[1].split('\n@<TRIPOS>BOND')[0].split("\n")
                    self.bonds = ''.join(lines).split('@<TRIPOS>BOND\n')[1].split("\n")
        
                    self.new_atoms = []
                    self.new_bonds = []
                    self.populate_mol(path=molecule_mol2, lines=lines)
                    
                    self.new_atoms.clear()
                    self.new_bonds.clear()
                    self.populate_mol(path=cell_mol2, lines=lines, seed=''.join(lines).split('@<TRIPOS>ATOM\n')[1].split("\n")[:len(self.atom_positions)])
            else:
                self.printToLog("# WARN - Compound ["+self.refcode+"] .mol2 output not found")
        except subprocess.CalledProcessError as e:
            self.printToLog("# WARN - Compound ["+self.refcode+"] Error creating _super.mol2")
            self.printToLog(str(e))
    
    
class featureExtractor:
    def printToLog(self,info):#Prints and logs in one, convention I personally like
        printToLog(self.log, info)    
    
    def getAttached(self,atom):
        atoms = list(map(lambda bond: bond.GetBeginAtom() if bond.GetEndAtomIdx() == atom.GetIdx() else bond.GetEndAtom(), atom.GetBonds()))
        return sorted(atoms, key=lambda x: -int(x.GetAtomicNum()))
        
    def parseCoord(self,coords):
        return [coords.x,coords.y,coords.z]
    
    def getOpposite(self,conf, atom, centre):
        return next(iter(sorted(self.getAttached(centre), key=lambda x: -int(GetAngleDeg(conf,x.GetIdx(),centre.GetIdx(),atom.GetIdx())))), None)

    def sigma(self, index):
        if any(line.startswith(str(self.atoms[index].strip().split()[-1])) for line in self.summary_atoms):
            for number, line in enumerate(self.summary_atoms, 0):
                if line.startswith(str(self.atoms[index].strip().split()[-1])):
                    self.data[f"{self.mol.GetAtoms()[index].GetSymbol()} [{index}] sigma"] = line.split("(")[1].split(")")[0]
                    self.data[f"{self.mol.GetAtoms()[index].GetSymbol()} [{index}] sigma_11"] = line.split("[sigma_11 ")[1].split("]")[0]
                    self.data[f"{self.mol.GetAtoms()[index].GetSymbol()} [{index}] sigma_22"] = line.split("[sigma_22 ")[1].split("]")[0]
                    self.data[f"{self.mol.GetAtoms()[index].GetSymbol()} [{index}] sigma_33"] = line.split("[sigma_33 ")[1].split("]")[0]
        else:
            self.printToLog(f"# WARN - No sigma values associated with atom [{self.mol.GetAtoms()[index].GetSymbol()}] [{index}] attatched to [{self.refcode}_{self.site_id}]")
            self.data[f"{self.mol.GetAtoms()[index].GetSymbol()} [{index}] sigma"] = "Not Found"
            self.data[f"{self.mol.GetAtoms()[index].GetSymbol()} [{index}] sigma_11"] = "Not Found"
            self.data[f"{self.mol.GetAtoms()[index].GetSymbol()} [{index}] sigma_22"] = "Not Found"
            self.data[f"{self.mol.GetAtoms()[index].GetSymbol()} [{index}] sigma_33"] = "Not Found"

    def __init__(self, log, directory, refcode, summary_atoms, site_id):
        self.log = log
        self.directory = directory
        self.refcode = refcode
        self.site_id = site_id
        self.summary_atoms = summary_atoms
        self.data = {}

        self.bond_cutoff = 2
        
    def extract(self):        
        with open(os.path.join(self.directory, self.refcode+"_cell.mol2")) as cell:
            self.atoms = ''.join(cell.readlines()).split('@<TRIPOS>ATOM\n')[1].split('\n@<TRIPOS>BOND')[0].split("\n")
            for W_number, W_line in enumerate(self.atoms, 0):
                if self.site_id in W_line:
                    W_id = W_number
                    break
                            
        self.mol = Chem.MolFromMol2File(os.path.join(self.directory, self.refcode+"_cell.mol2"), sanitize=False, removeHs=False)
        if self.mol == None:
            return None
        
        ComputeCanonicalTransform(self.mol.GetConformer(),center=self.mol.GetConformer().GetAtomPosition(W_id))
        conf = self.mol.GetConformer()
        W = self.mol.GetAtoms()[W_id]
        self.printToLog(f"# INFO - Compound [{self.refcode}_{self.site_id}] Processing site [{W.GetSymbol()}] [{W_id}]")


        
        W_attached = self.getAttached(W)
        symbols = list(map(lambda atom: atom.GetSymbol(), W_attached))                                
        expected = ["P", "C", "C", "C", "C", "C"]
        
        self.printToLog(f"# INFO - Compound [{self.refcode}_{self.site_id}] Atom [{W.GetSymbol()}] [{W_id}] has the following neighbours [{symbols}]")
        if len(W_attached) == 6 and symbols == expected:   
            self.data[f"Site ID"] = f"{self.refcode}_{self.site_id}"
            
            self.sigma(W.GetIdx())

            P = list(filter(lambda x: x.GetSymbol() == "P", W_attached))[0]
            self.sigma(P.GetIdx())

            W_site_coords = [self.parseCoord(conf.GetAtomPosition(W_id))]
            W_site_coords.extend(list(map(lambda atom: self.parseCoord(conf.GetAtomPosition(atom.GetIdx())), W_attached)))
            
            dist = oc.CalcDistortion(W_site_coords)
            self.data[f"W Distortion zeta"] = dist.zeta
            self.data[f"W Distortion delta"] = dist.delta
            self.data[f"W Distortion sigma"] = dist.sigma
            self.data[f"W Distortion theta"] = dist.theta

            carbons = sorted(W_attached[1:].copy(), key=lambda x: x.GetIdx() == self.getOpposite(conf, P, W).GetIdx())
            processed = []
            for C in carbons:                                    
                opposite = self.getOpposite(conf, C, W)
                if processed.__contains__(C.GetIdx()) or processed.__contains__(opposite.GetIdx()):
                    continue

                processed.append(C.GetIdx())
                processed.append(opposite.GetIdx())

                # Angle between opposite sides of W
                self.data[f"{C.GetSymbol()} {W.GetSymbol()} {opposite.GetSymbol()} [{C.GetIdx()}] [{W.GetIdx()}] [{opposite.GetIdx()}]"] = GetAngleDeg(conf,C.GetIdx(),W.GetIdx(),opposite.GetIdx())

                # Distance to W
                self.data[f"{C.GetSymbol()} {W.GetSymbol()} [{C.GetIdx()}] [{W.GetIdx()}]"] = GetBondLength(conf,C.GetIdx(),W.GetIdx())
                
                # Bonds to atom bonded to W
                for O in self.getAttached(C)[1:]:
                    self.data[f"{C.GetSymbol()} {O.GetSymbol()} [{C.GetIdx()}] [{O.GetIdx()}]"] = GetBondLength(conf,C.GetIdx(),O.GetIdx())
                    self.data[f"{W.GetSymbol()} {C.GetSymbol()} {O.GetSymbol()} [{W.GetIdx()}] [{C.GetIdx()}] [{O.GetIdx()}]"] = GetAngleDeg(conf,W.GetIdx(),C.GetIdx(),O.GetIdx())

                self.sigma(C.GetIdx())
                
   
                
                if not opposite.GetIdx() == P.GetIdx():
                    self.data[f"{opposite.GetSymbol()} {W.GetSymbol()} [{opposite.GetIdx()}] [{W.GetIdx()}]"] = GetBondLength(conf,opposite.GetIdx(),W.GetIdx())
                    for O in self.getAttached(opposite)[1:]:
                        self.data[f"{opposite.GetSymbol()} {O.GetSymbol()} [{opposite.GetIdx()}] [{O.GetIdx()}]"] = GetBondLength(conf,opposite.GetIdx(),O.GetIdx())
                        self.data[f"{W.GetSymbol()} {opposite.GetSymbol()} {O.GetSymbol()} [{W.GetIdx()}] [{opposite.GetIdx()}] [{O.GetIdx()}]"] = GetAngleDeg(conf,W.GetIdx(),opposite.GetIdx(),O.GetIdx())
                    self.sigma(opposite.GetIdx())
                
   
            
            P_attached = self.getAttached(P)
            symbols = list(map(lambda atom: atom.GetSymbol(), P_attached))                                
            self.printToLog(f"# INFO - Compound [{self.refcode}_{self.site_id}] Atom [{W.GetSymbol()}] [{W_id}] attached P has the following neighbours [{symbols}]")

            if len(P_attached) <= self.bond_cutoff or len(P_attached) > 4:
                self.printToLog(f"# WARN - Compound [{self.refcode}_{self.site_id}] Atom [{W.GetSymbol()}] [{W_id}] attached P has [{len(P_attached)}] P bonds")
                return None
                
            P_site_dist = list(map(lambda atom: GetBondLength(conf,P.GetIdx(),atom.GetIdx()), P_attached))
            P_site_dist_mean = sum(dist for dist in P_site_dist) / len(P_site_dist)

            self.data[f"P Distortion zeta"] = sum(abs(dist - P_site_dist_mean) for dist in P_site_dist)
            self.data[f"P Distortion delta"] = sum(pow((dist - P_site_dist_mean) / P_site_dist_mean, 2) for dist in P_site_dist) / len(P_site_dist)
            P_site_angle = []
            for i in range(len(P_attached)):
                for j in range(i + 1, len(P_attached)):
                    P_site_angle.append(GetAngleDeg(conf,P_attached[i].GetIdx(),P.GetIdx(),P_attached[j].GetIdx()))
            self.data[f"P Distortion sigma"] = sum(abs(109.5 - angle) for angle in P_site_angle)

            self.data[f"Bonds {P.GetSymbol()} {P.GetIdx()}"] = len(P_attached)
            self.data[f"{P.GetSymbol()} {W.GetSymbol()} [{P.GetIdx()}] [{W.GetIdx()}]"] = GetBondLength(conf,P.GetIdx(),W.GetIdx())

            # Bonds to P bonded to W
            for atom in P_attached:
                if not atom.GetIdx() == W.GetIdx():
                    self.data[f"Atomic Number {atom.GetSymbol()} {atom.GetIdx()}"] = atom.GetAtomicNum()
                    self.data[f"Bonds {atom.GetSymbol()} {atom.GetIdx()}"] = int(len(self.getAttached(atom)))
                    self.data[f"{W.GetSymbol()} {P.GetSymbol()} {atom.GetSymbol()} [{W.GetIdx()}] [{P.GetIdx()}] [{atom.GetIdx()}]"] = GetAngleDeg(conf,W.GetIdx(),P.GetIdx(),atom.GetIdx())
                    self.data[f"{P.GetSymbol()} {atom.GetSymbol()} [{P.GetIdx()}] [{atom.GetIdx()}]"] = GetBondLength(conf,P.GetIdx(),atom.GetIdx())

            return self.data
        else:
            self.printToLog(f"# WARN - Compound [{self.refcode}_{self.site_id}] Atom [{W.GetSymbol()}] [{W_id}] is invalid")
            return None