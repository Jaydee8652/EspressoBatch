import re 
import os
import sys
import time
import datetime
import pandas as pd
import numpy as np
import datetime
import time

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.svm import SVR

from utils.generic_utils import printToLog as pl, createDirectory as cd, cellVolume
from utils.params import *

from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor, rdFMCS, rdDistGeom
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeProperties, MMFFGetMoleculeForceField
from rdkit.Chem.rdMolTransforms import *

import optuna
import octadist as oc
from sklearn.multioutput import MultiOutputRegressor

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)
    
def getAttached(atom):
    atoms = list(map(lambda bond: bond.GetBeginAtom() if bond.GetEndAtomIdx() == atom.GetIdx() else bond.GetEndAtom(), atom.GetBonds()))
    return sorted(atoms, key=lambda x: -int(x.GetAtomicNum()))
def parseCoord(coords):
    return [coords.x,coords.y,coords.z]
    
#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

printToLog("# INFO - Enter integer to choose process to perform.")
options = {
    "1": "Recreate sheet",
    "2": "Run GradientBoostingRegressor - Multiple Output",
    "3": "Run GradientBoostingRegressor - Single Output"
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

forest_data = os.path.join(homeDirectory, "_forest_data.csv")
if choices.__contains__("1"):
    if os.path.isfile(forest_data):
        os.remove(forest_data)
    printToLog("# INFO - Creating new sheet ["+ forest_data + "]")
    with open(forest_data, 'a') as file:
        print("REFCODE,W sigma,W sigma_11,W sigma_22,W sigma_33,P sigma,P sigma_11,P sigma_22,P sigma_33,Distortion zeta,Distortion delta,Distortion sigma,Distortion theta,P bonds,P-W dist,?1 mass,?1 bonds,W-P-?1 angle,P-?1 dist,?2 mass,?2 bonds,W-P-?2 angle,P-?2 dist,?3 mass,?3 bonds,W-P-?3 angle,P-?3 dist", file=file)

    
    Output_Files = os.path.join(homeDirectory, "Output_Files")
    createDirectory(Output_Files, "# WARN - No directory found for input files.", True)
    directories = sorted([directory for directory in os.listdir(Output_Files) if os.path.isdir(os.path.join(Output_Files, directory))])
    
    if len(directories) == 0:
        printToLog("# WARN - No directories found in ["+ Output_Files + "]")
        quit()
    else:
        print(directories)
        for refcode in directories:
            printToLog(f"# INFO - Processing compound [{refcode}]")
            
            refcodeDirectory = os.path.join(Output_Files, refcode)
            summary = os.path.join(refcodeDirectory, refcode+"_summary.txt")
            
            if os.path.exists(summary):
                with open(summary) as file:
                    summary_lines = file.readlines()
                    
                    join = ''.join(summary_lines)
                    if "BATCH_done = True" in join and "PWSCF_done = True" in join and "PWSCF_finalEnergy" in join and "GIPAW_done = True" in join:        
                        summary_atoms = join.split('#BEGIN_ATOMIC_POSITIONS\n')[1].split('\n#END_ATOMIC_POSITIONS')[0].split("\n")
                        for summary_number, summary_line in enumerate(summary_atoms, 0):
                            if summary_line.split()[1] == "W":
                                site_id = summary_line.split()[0]
                                
                                with open(os.path.join(refcodeDirectory, refcode+"_cell.mol2")) as cell:
                                    atoms = ''.join(cell.readlines()).split('@<TRIPOS>ATOM\n')[1].split('\n@<TRIPOS>BOND')[0].split("\n")
                                    for W_number, W_line in enumerate(atoms, 0):
                                        if site_id in W_line:
                                            W_id = W_number
                                            break
                                                    
                                mol = Chem.MolFromMol2File(os.path.join(refcodeDirectory, refcode+'_cell.mol2'), sanitize=False, removeHs=False)
                                if mol == None:
                                    continue
                                
                                ComputeCanonicalTransform(mol.GetConformer(),center=mol.GetConformer().GetAtomPosition(W_id))
                                conf = mol.GetConformer()
                                W = mol.GetAtoms()[W_id]
                                printToLog(f"# INFO - Compound [{refcode}] Processing atom W [#{W_id} {W.GetSymbol()}]")
    
                                data = {}
                                P_BOND_COUNT = 0
    
                                local = getAttached(W)
                                symbols = list(map(lambda atom: atom.GetSymbol(), local))                                
                                expected = ["P", "C", "C", "C", "C", "C"]
                                
                                printToLog(f"# WARN - Compound [{refcode}] Atom W [#{W_id}] has the following neighbours [{symbols}]")
                                if len(local) == 6 and symbols == expected:   
                                    data[f"Site ID"] = f"{refcode}_{site_id}"
                                    data[f"W sigma"] = summary_line.split("(")[1].split(")")[0]
                                    data[f"W sigma_11"] = summary_line.split("[sigma_11 ")[1].split("]")[0]
                                    data[f"W sigma_22"] = summary_line.split("[sigma_22 ")[1].split("]")[0]
                                    data[f"W sigma_33"] = summary_line.split("[sigma_33 ")[1].split("]")[0]

                                    P = list(filter(lambda x: x.GetSymbol() == "P", local))[0]
                                    for P_number, P_line in enumerate(summary_atoms, 0):
                                        if P_line.startswith(str(atoms[P.GetIdx()].strip().split()[-1])):
                                            data[f"P sigma"] = P_line.split("(")[1].split(")")[0]
                                            data[f"P sigma_11"] = P_line.split("[sigma_11 ")[1].split("]")[0]
                                            data[f"P sigma_22"] = P_line.split("[sigma_22 ")[1].split("]")[0]
                                            data[f"P sigma_33"] = P_line.split("[sigma_33 ")[1].split("]")[0]
                                            break

                                    coord = [parseCoord(conf.GetAtomPosition(W_id))]
                                    coord.extend(list(map(lambda atom: parseCoord(conf.GetAtomPosition(atom.GetIdx())), local)))
                                    
                                    dist = oc.CalcDistortion(coord)
                                    zeta = dist.zeta
                                    delta = dist.delta
                                    sigma = dist.sigma
                                    theta = dist.theta

                                    data[f"Distortion zeta"] = dist.zeta
                                    data[f"Distortion delta"] = dist.delta
                                    data[f"Distortion sigma"] = dist.sigma
                                    data[f"Distortion theta"] = dist.theta
                                                                                                        
                                    P_attached = getAttached(P)
                                    data[f"Bonds {P.GetSymbol()} {P.GetIdx()}"] = len(P_attached)
                                    
                                    # Distance to W
                                    data[f"{P.GetSymbol()} {W.GetSymbol()} [{P.GetIdx()}] [{W.GetIdx()}]"] = GetBondLength(conf,P.GetIdx(),W.GetIdx())

                                    # Bonds to P bonded to W
                                    for atom in P_attached:
                                        if not atom.GetIdx() == W.GetIdx():
                                            data[f"Mass {atom.GetSymbol()} {atom.GetIdx()}"] = atom.GetAtomicNum()
                                            data[f"Bonds {atom.GetSymbol()} {atom.GetIdx()}"] = int(len(getAttached(atom)))
                                            data[f"{W.GetSymbol()} {P.GetSymbol()} {atom.GetSymbol()} [{W.GetIdx()}] [{P.GetIdx()}] [{atom.GetIdx()}]"] = GetAngleDeg(conf,W.GetIdx(),P.GetIdx(),atom.GetIdx())
                                            data[f"{P.GetSymbol()} {atom.GetSymbol()} [{P.GetIdx()}] [{atom.GetIdx()}]"] = GetBondLength(conf,P.GetIdx(),atom.GetIdx())

                                                
                                    values = list(data.values())
                                    df = pd.read_csv(forest_data)
                                    while len(values) < len(df.columns):
                                        values.append("0")
                                    df.loc[len(df)] = values
                                    
                                    df.to_csv(forest_data, index=False)
                                else:
                                    printToLog(f"# WARN - Compound [{refcode}] Atom W [{W_id}] invalid")
    
                    else:
                        printToLog(f"# WARN - Compound [{refcode}] Not complete")
            else:
                printToLog(f"# WARN - Compound [{refcode}] No summary file found")                            

if not os.path.isfile(forest_data):
    printToLog("# WARN - No .csv found.")
    quit()

df = pd.read_csv(forest_data)
if choices.__contains__("2"):
    df = df.drop(['W sigma'], axis=1)
    X = df.drop(['REFCODE','W sigma_11', 'W sigma_22', 'W sigma_33'], axis=1)
    y = df[['REFCODE','W sigma_11','W sigma_22','W sigma_33']]
elif choices.__contains__("3"):
    df = df.drop(['W sigma_11', 'W sigma_22', 'W sigma_33'], axis=1)
    X = df.drop(['REFCODE', 'W sigma'], axis=1)
    y = df[['REFCODE','W sigma']]

def randomForest(loss='squared_error', learning_rate=0.1, n_estimators=100, subsample=1.0, min_samples_split=2, min_samples_leaf=1, min_weight_fraction_leaf=0.0, max_depth=3, min_impurity_decrease=0.0, init=None, random_state=42, max_features=None, alpha=0.9, verbose=0, max_leaf_nodes=None, warm_start=False, validation_fraction=0.1, n_iter_no_change=None, tol=0.0001, ccp_alpha=0.0):
    forest_paramaters = locals()
    printToLog(f"# INFO - Running Random Forest with paramaters [{forest_paramaters}]")
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=random_state)
    y_train = y_train.drop('REFCODE', axis=1)
    
    regressor = GradientBoostingRegressor(**forest_paramaters)
    if choices.__contains__("2"):
        regressor = MultiOutputRegressor(regressor)
 
    regressor.fit(X_train, y_train)
    
    y_pred = regressor.predict(X_test)
    
    mse = mean_squared_error(y_test.drop('REFCODE', axis=1), y_pred)
    r2 = r2_score(y_test.drop('REFCODE', axis=1), y_pred)

    predicton_data = os.path.join(prediction, f"_predicton_data_{random_state}.csv")
    if os.path.isfile(predicton_data):
        os.remove(predicton_data)
    printToLog("# INFO - Creating new sheet ["+ predicton_data + "]")

    with open(predicton_data, 'a') as file:
        if choices.__contains__("2"):
            print("REFCODE,ACTUAL_sigma_11,PREDICTED_sigma_11,DIFFERENCE_sigma_11,ACTUAL_sigma_22,PREDICTED_sigma_22,DIFFERENCE_sigma_22,ACTUAL_sigma_33,PREDICTED_sigma_33,DIFFERENCE_sigma_33",file=file)
            for i in range(len(y_pred)):
                refcode = str(y_test.iloc[i]['REFCODE'])
                pred_df = pd.DataFrame(y_pred)
    
                actual_sigma_11 = float(y_test.iloc[i]['W sigma_11'])
                pred_sigma_11 = pred_df[0][i]
                difference_sigma_11 = abs(actual_sigma_11-pred_sigma_11)            
    
                actual_sigma_22 = float(y_test.iloc[i]['W sigma_22'])
                pred_sigma_22 = pred_df[1][i]
                difference_sigma_22 = abs(actual_sigma_22-pred_sigma_22)            
    
                actual_sigma_33 = float(y_test.iloc[i]['W sigma_33'])
                pred_sigma_33 = pred_df[2][i]
                difference_sigma_33 = abs(actual_sigma_33-pred_sigma_33)            
                
                print(f"{refcode},{actual_sigma_11:.2f},{pred_sigma_11:.2f},{difference_sigma_11:.2f},{actual_sigma_22:.2f},{pred_sigma_22:.2f},{difference_sigma_22:.2f},{actual_sigma_33:.2f},{pred_sigma_33:.2f},{difference_sigma_33:.2f}",file=file)
        elif choices.__contains__("3"):
            print("REFCODE,ACTUAL,PREDICTED,DIFFERENCE",file=file)
            for i in range(len(y_pred)):
                refcode = str(y_test.iloc[i]['REFCODE'])
                actual = float(y_test.iloc[i]['W sigma'])
                predicted = float(y_pred[i])
                difference = abs(actual-predicted)            
                print(f"{refcode},{actual:.2f},{predicted:.2f},{difference}",file=file)
            
        print(f"Paramaters,{forest_paramaters}",file=file)
        print(f"Mean Squared Error,{mse:.2f}",file=file)
        print(f"R-squared Score,{r2:.2f}",file=file)
    
    with open(summary, 'a') as file:
        print(f"{mse:.2f},{r2:.2f},{forest_paramaters}", file=file)
    
    printToLog(f"# INFO - Mean Squared Error: {mse:.2f}")
    printToLog(f"# INFO - R-squared Score: {r2:.2f}")

def forest_objective(trial, random_state=42):
    param_space = {
        'n_estimators': trial.suggest_int('n_estimators', 50, 500, step=25),
        'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.2),
        'max_depth': trial.suggest_int('max_depth', 1, 6),
    }
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=random_state)
    y_train = y_train.drop('REFCODE', axis=1)

    regressor = GradientBoostingRegressor(**param_space, random_state=random_state)
    if choices.__contains__("2"):
        regressor = MultiOutputRegressor(regressor)
 
    regressor.fit(X_train, y_train)
    y_pred = regressor.predict(X_test)
    return r2_score(y_test.drop('REFCODE', axis=1), y_pred) 

if choices.__contains__("2") or choices.__contains__("3"):
    out_type = "single"
    if choices.__contains__("2"):
        out_type = "multi"
        
    time = str(datetime.datetime.now().strftime("[%Y-%m-%d_%H-%M-%S]"))
    prediction = os.path.join(homeDirectory, f"prediction_{out_type}_{time}")
    createDirectory(prediction, "# WARN - No directory found for prediction files.", False)
    summary = os.path.join(prediction, "prediction_summary.csv")
    if os.path.isfile(summary):
        os.remove(summary)
    with open(summary, 'a') as file:
        print(f"MSE,R-SQUARED,PARAMS", file=file)
        
    for random_state in range(100):
        study = optuna.create_study(direction='maximize')
        study.optimize(lambda trial: forest_objective(trial, random_state=random_state), n_trials=50)
        printToLog(f"# INFO - Best paramaters determined to be [{study.best_params}]")
        randomForest(**study.best_params, random_state=random_state)
