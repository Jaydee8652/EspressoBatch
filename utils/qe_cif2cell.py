# qe_cif2cell.py - Jacob Duddridge

# Generates quantumespresso and gipaw input files from a .cif
# Runs test calculation and updates resource requests accordingly
# Atom fixing is included

# All processes are reported to qe_cif2cell.log for debugging

#Imports
import __main__ as main
import os
import subprocess
import sys
import shutil
import re 
import math
import datetime
import time
import pandas as pd

#Jank thing to fix the path. Very annoying artefact of running python scripts by absolute path.
sys.path[0] = sys.path[0][:-6] + sys.path[0][-6:].replace("/utils", "")

from utils.generic_utils import printToLog as pl, createDirectory as cd, cellVolume
from utils.params import *

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)
    
#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")

refcodeDirectory = os.getcwd() #Directory where we are
input_path = os.path.split(refcodeDirectory)[0]
set_directory = os.path.split(input_path)[0]
homeDirectory = os.path.split(set_directory)[0]

refcode = os.path.basename(refcodeDirectory)

if len(sys.argv) > 1:
    printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

local_params = os.path.join(set_directory, "local_params.csv")
df = pd.read_csv(local_params, encoding="utf-8-sig")

#params
set_id = str(df['set_id'].iloc[0])
ecutwfc = float(df['ecutwfc'].iloc[0])
ecutrho_factor = float(df['ecutrho_factor'].iloc[0])
conv_thr = str(df['conv_thr'].iloc[0])
calculation = str(df['calculation'].iloc[0])
q_gipaw = float(df['q_gipaw'].iloc[0])
volume_cap = float(df['volume_cap'].iloc[0])
atoms_to_optimise = str(df['atoms_to_optimise'].iloc[0])

#Make sure there is a directory to process
cifs_path = os.path.join(homeDirectory, "cifs")
createDirectory(cifs_path, "# WARN - No directory found for .cifs to process.", False)

validated = os.path.join(cifs_path, "validated")
createDirectory(validated, "# WARN - No directory found for .cifs to process. Place .cif files in or replace the newly created directory at ["+validated+"]", True)

cif = os.path.join(validated, refcode+".cif")
qe_in = os.path.join(refcodeDirectory, refcode+".in")
gipaw_in = os.path.join(refcodeDirectory, f"gipaw.{refcode}.in")
test_out = os.path.join(refcodeDirectory, f"test_{refcode}.out")

QE_SUB = os.path.join(refcodeDirectory, "QE_SUB")    
post = os.path.join(os.path.join(homeDirectory,"utils"), "post_processing.py")
incomplete = os.path.join(refcodeDirectory, "INCOMPLETE.txt")



if len(sys.argv) > 1:
    #Make sure there is a directory for PSEUDOS
    PSEUDOS = os.path.join(os.path.join(os.path.join(homeDirectory, "utils"), "data"), "PSEUDOS/")
    createDirectory(PSEUDOS, "# WARN - No directory found for PSEUDOS. Place .UPF files in or replace the newly created directory at", True)

    #Get .UPFs from directory
    pseudos = [file for file in os.listdir(PSEUDOS) if file.endswith('.UPF') and os.path.isfile(os.path.join(PSEUDOS, file))]
    
    if len(pseudos) == 0:#Make sure there are .UPFs in the directory
        printToLog("# WARN - No .UPF files found. Place .UPF files in ["+ PSEUDOS + "]")
        quit()
    else:
        printToLog("# INFO - " + str(len(pseudos)) + " .UPF files found at ["+ PSEUDOS + "]")
        psuedo_elements = [file.split(".")[0] for file in pseudos]
        printToLog("# INFO - Following elements accounted for ["+str(list(set(psuedo_elements)))+"]")

    # Input sanitisation - user input of atoms to freeze
    regex = re.compile('[^a-zA-Z ]')
    atoms_to_optimise = regex.sub('', atoms_to_optimise).strip().split(" ")
    atoms_to_optimise = list(set([atom.lower().capitalize() for atom in atoms_to_optimise]))

    if atoms_to_optimise.__contains__("All"): # If there is "All", remove everything else. "H C All" -> "All"
        atoms_to_optimise.clear()
        atoms_to_optimise.append("All")
    if atoms_to_optimise.__contains__("None"): # If there is "None", remove everything. "H C None" -> ""
        atoms_to_optimise.clear()

    invalidInputs = []
    for atom in atoms_to_optimise:    
        if not psuedo_elements.__contains__(atom) and not atom == "All":
            invalidInputs.append(atom)
    
    # If we don't have the pseudo for the atom, scrap it, we couldn't run it anyway
    if len(invalidInputs) > 0:
        printToLog("# WARN - The following atom types ["+str(list(set(invalidInputs)))+"] are not accounted for by the available .UPF files and have been removed")
    atoms_to_optimise = [atom for atom in atoms_to_optimise if atom not in invalidInputs]
    printToLog("# INFO - Following atom types selected to be optimised ["+str(atoms_to_optimise)+"]")

    printToLog("# INFO - Compound [" + refcode + "] Processing .cif file")
    shutil.copyfile(cif, os.path.join(refcodeDirectory, refcode+".cif"))

    try:
        subprocess.run(f"cif2cell -f {cif} -p quantum-espresso --setup-all -o {qe_in}",shell=True)
        printToLog("# INFO - Compound [" + refcode + "] Sucessfully ran cif2cell")
    except subprocess.CalledProcessError as e:
        printToLog("# WARN - Compound [" + refcode + "] Error running cif2cell")
        printToLog(str(e))
        quit()
    
    if os.path.exists(qe_in):
        with open(qe_in) as file:
            lines = file.readlines()
            offset = 0
            for number, line in enumerate(lines.copy(), 0): 
                number += offset
    
                line_atom = line.lstrip().split(" ")[0]
                # Freeze atoms
                if not atoms_to_optimise.__contains__(line_atom) and not atoms_to_optimise.__contains__("All") and psuedo_elements.__contains__(line_atom) and len(line.split()) > 3:
                    lines[number] = lines[number].rstrip() + " 0 0 0\n"            
                
                if "&SYSTEM" in line:
                    content = f"""
&CONTROL
  calculation = '{calculation}'
  prefix = '{refcode}'
  outdir = '{refcodeDirectory}/'
  pseudo_dir = '{PSEUDOS}'
  nstep = 0
/

"""
                    lines.insert(number, content.lstrip("\n"))
                    offset += 1
                if "ntyp" in line:
                    del lines[number+1]
                    content = f"""
  ecutwfc = {ecutwfc}
  ecutrho = {ecutwfc * ecutrho_factor}
  vdw_corr = 'grimme-d3'
  dftd3_version = 6
/

&ELECTRONS
  conv_thr = {conv_thr}
/

&IONS
/

"""
                    lines.insert(number+1, content.lstrip("\n"))
                if "ATOMIC_SPECIES" in line:
                    lines.insert(number,"\n")
                    offset += 1
                if "ATOMIC_POSITIONS" in line:
                    lines.insert(number,"\n")
                    offset += 1
                if "_PSEUDO" in line:
                    line_atom = line.split()[0]
                    if line_atom+"_PSEUDO" in line:                # Add pseuds
                        if any(line_atom+'.' in pseudo and 'kjpaw' in pseudo for pseudo in pseudos):
                            for pseudo in pseudos:
                                if line_atom+'.' in pseudo and 'kjpaw' in pseudo:
                                    lines[number]=line.replace(line_atom+"_PSEUDO",pseudo)
                        else:
                            printToLog("# WARN - No .UPF file found for atom ["+str(line_atom.lstrip())+"], needed by compound ["+refcode+"]")
                            shutil.rmtree(refcodeDirectory)
                            quit() # Deliberate crash if a pseudopotential is not available for all atom types required
                            
        with open(qe_in,"w") as file:
            for line in lines:
                file.write(line)
        printToLog("# INFO - Compound [" + refcode + "] Created .in file at ["+qe_in+"]")
    
        # Create the .gipaw input file
        with open(gipaw_in, "w") as file:
            content = f"""
&inputgipaw
  job = 'nmr'
  prefix = '{refcode}'
  tmp_dir = '{refcodeDirectory}/'
  diagonalization = 'david'
  verbosity = 'high'
  q_gipaw = {q_gipaw}
  spline_ps = .true.
  use_nmr_macroscopic_shape = .true.
/        
"""
            print(content.lstrip(), file=file)
        printToLog("# INFO - Compound [" + refcode + "] Created gipaw .in file at ["+gipaw_in+"]")
        printToLog("# INFO - Compound [" + refcode + "] Running test calculation")

else:
    # Process test result
    if os.path.exists(test_out):
        printToLog("# INFO - Compound [" + refcode + "] Processing test output file")
        
        kPoints = 0
        dynamicalRAM = 0
        
        with open(test_out, "r") as file:
            for line in file:                
                if len(line) == 0 or line.startswith("#"):
                    continue # skip blank lines and comments
                if "number of k points" in line:
                    kPoints = math.ceil(float(re.sub("[^0-9]", "", line)))
                    printToLog("# INFO - Compound [" + refcode + "] Determined to have ["+str(kPoints)+"] k points")    
                if "Estimated max dynamical RAM per process" in line:
                    dynamicalRAM = (math.ceil(float(re.sub("[^0-9.]", "", line))) * 3)
                    printToLog("# INFO - Compound [" + refcode + "] Determined to use ["+str(dynamicalRAM)+"G] max dynamical RAM per process")    
                    
        if dynamicalRAM > (3.7 * kPoints):
            mult = math.ceil(dynamicalRAM / (3.7 * kPoints))
            kPoints *= mult
            printToLog("# INFO - Compound [" + refcode + "] Has high RAM usage. Multiplied number of tasks per node by ["+str(mult)+"]")
        
        if kPoints > param_cores:
            printToLog("# WARN - Compound [" + refcode + "] Number of tasks per node ["+str(kPoints)+"] above cap of ["+str(param_cores)+"]")   
            kPoints = param_cores
        if dynamicalRAM > param_memory:
            printToLog("# WARN - Compound [" + refcode + "] Projected RAM usage ["+str(dynamicalRAM)+"G] above cap of ["+str(param_memory)+"G]")   
            dynamicalRAM = param_memory

        with open(cif, "r") as file:
            lines = file.readlines()
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
        
        days = "00"
        if volume > volume_cap:
            printToLog(f"# INFO - Compound [{refcode}] Volume greater than [{volume_cap}], 2 extra days allocated.")
            days = "02"
        
        printToLog(f"# INFO - Compound [{refcode}] Number of tasks per node set to [{kPoints}]")
        printToLog(f"# INFO - Compound [{refcode}] Memory request set to [{dynamicalRAM}G]")

        # Create the QE_SUB submission script
        with open(QE_SUB, "w") as file:
            content = f"""
#!/bin/bash
    
#SBATCH --job-name={set_id}_{refcode}_SUB
#SBATCH --mail-type={param_slurmVerbosity}
#SBATCH --mail-user={param_email}
#SBATCH --account={param_account}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={kPoints}
#SBATCH --cpus-per-task=1
#SBATCH --mem={dynamicalRAM}G
#SBATCH --time={days}-23:59

module load {param_modules}

srun --cpus-per-task=$SLURM_CPUS_PER_TASK pw.x < {refcode}.in > {refcode}.out
srun --cpus-per-task=$SLURM_CPUS_PER_TASK gipaw.x < gipaw.{refcode}.in > gipaw.{refcode}.out
srun --cpus-per-task=1 --ntasks=1 python3 {post}
"""
            print(content.lstrip("\n"),file=file)
        printToLog("# INFO - Compound [" + refcode + "] Created QE_SUB file at ["+QE_SUB+"]")
    
        printToLog("# INFO - Compound [" + refcode + "] Updating .in file")
        with open(qe_in) as file:
            lines = file.read().splitlines()
        with open(qe_in, "w") as file:
            for line in lines:
                if "nstep" in line:
                    print(re.sub("  nstep = 0", "  !nstep = 0", line), file=file)
                    printToLog("# INFO - Compound [" + refcode + "] Ready to run real calculation!")
                    os.remove(incomplete)
                else:
                    print(line, file=file)
    else:
        printToLog("# WARN - Compound [" + refcode + "] No test output found")
        quit()