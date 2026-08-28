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

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from matplotlib import colormaps
from PIL import Image

import joblib

from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.svm import SVR
from sklearn.inspection import permutation_importance

from utils.generic_utils import printToLog as pl, createDirectory as cd, cellVolume
from utils.params import *

from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor, rdFMCS, rdDistGeom
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeProperties, MMFFGetMoleculeForceField
from rdkit.Chem.rdMolTransforms import *

import optuna
import octadist as oc
from sklearn.multioutput import MultiOutputRegressor

n_trials = 100 # 50
n_states = 100 # 100
n_repeats = 150 # 150
cmap = plt.cm.viridis

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
def getOpposite(conf, atom, centre):
    return next(iter(sorted(getAttached(centre), key=lambda x: -int(GetAngleDeg(conf,x.GetIdx(),centre.GetIdx(),atom.GetIdx())))), None)

def plotPrediction(random_state):
    title = f"_prediction_data_{random_state}"
    printToLog(f"# INFO - Plotting [{title}]")

    plt.clf()
    prediction_data = os.path.join(prediction, f"{title}.csv")
    df = pd.read_csv(prediction_data)
    df.set_index("REFCODE", inplace = True)
    df = df.drop(["Parameters","Mean Squared Error","R-squared Score"])
    df = df.astype(float)
    headers = list(df.columns.values)
    
    
    Actual = df.iloc[:, 0]    
    Predicted = df.iloc[:, 1]
    
    plt.scatter(Actual, Predicted, s=1)
    plt.title(f"{title}")
    
    axes = plt.gca()
    axes.set_xlabel(headers[0])
    axes.set_ylabel(headers[1])
    
    xmin, xmax, ymin, ymax = plt.axis()
    plt.axline((xmin, ymin), slope=1, linewidth=1, color='red')
    plt.savefig(f'{plots}/{title}.png',bbox_inches='tight')
    
def plotPredictionSummary():
    title = f"prediction_summary"
    printToLog(f"# INFO - Plotting [{title}]")
    
    plt.clf()
    df = pd.read_csv(prediction_summary)
    df = df.drop(df.tail(1).index)
    df = df.drop(['Parameters'], axis=1)
    df = df.astype(float)
    headers = list(df.columns.values)
    
    x = df.iloc[:, 0]    
    y = df.iloc[:, 1]

    plt.scatter(x, y, s=1)
    plt.plot(x, [np.mean(y)]*len(x), linewidth=1, linestyle='dashed', color='red')
    plt.title(f"{title}")
    
    axes = plt.gca()

    axes.set_xlabel(headers[0])
    axes.set_ylabel(headers[1])
    
    plt.savefig(f'{plots}/{title}.png',bbox_inches='tight')

def errorHeatmap(x,y,xerr):
    plt.scatter(x,y,s=0.1,c=x,cmap=cmap)
    
    norm = matplotlib.colors.Normalize(vmin=0, vmax=0.5, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    error_color = np.array([(mapper.to_rgba(v)) for v in x])
    
    #loop over each data point to plot
    for x, y, e, color in zip(x, y, xerr, error_color):
        plt.plot(x, y, 'o',markersize=1, color=color)
        plt.errorbar(x,y,xerr=e,fmt='none',color=color,capsize=3,elinewidth=0.5,capthick=0.5)


def plotImportance(random_state):
    title = f"_importance_data_{random_state}"
    printToLog(f"# INFO - Plotting [{title}]")

    plt.clf()
    x = pd.read_csv(mean_importance_summary, index_col=["Random state"]).iloc[random_state]
    y = list(pd.read_csv(mean_importance_summary, index_col=["Random state"]).columns.values)
    xerr = pd.read_csv(deviation_importance_summary, index_col=["Random state"]).iloc[random_state]

    errorHeatmap(x,y,xerr)
    
    plt.title(f"{title}")
    
    axes = plt.gca()
    axes.set_xlim(0, 1)

    axes.set_xlabel("Mean Importance")
    axes.set_ylabel("Feature Name")

    plt.savefig(f'{plots}/{title}.png',bbox_inches='tight')

def plotImportanceSummary():
    title = f"importance_summary"
    printToLog(f"# INFO - Plotting [{title}]")

    plt.clf()

    mean = pd.read_csv(mean_importance_summary, index_col=["Random state"]).tail(1)
    mean.insert(loc=0, column="Index", value=["Mean"])
    mean = mean.set_index(["Index"])
    
    deviation = pd.read_csv(deviation_importance_summary, index_col=["Random state"]).tail(1)
    deviation.insert(loc=0, column="Index", value=["STD"])
    deviation = deviation.set_index(["Index"])

    df = pd.concat([mean,deviation])
    df = df.sort_values(by=["Mean"],axis=1,ascending=True)

    x = df.iloc[0, :]
    y = list(df.columns.values)
    xerr = df.iloc[1, :]

    errorHeatmap(x,y,xerr)
    
    plt.title(f"{title}")
    
    axes = plt.gca()    
    axes.set_xlabel("Mean Importance")
    axes.set_ylabel("Feature Name")

    plt.savefig(f'{plots}/{title}.png',bbox_inches='tight')

def plotCorrelation():
    title = f"feature_correlation"

    plt.matshow(X.corr(), cmap=cmap)
    plt.xticks(ticks = range(0,len(X.columns)), labels = list(X.columns.values), rotation=90)
    plt.yticks(ticks = range(0,len(X.columns)), labels = list(X.columns.values), rotation=0)
    plt.title(title)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14)
    plt.savefig(f'{plots}/{title}.png',bbox_inches='tight')
    

def createGrid(name):
    printToLog(f"# INFO - Creating grid of {name}")

    images = [Image.open(f"{plots}/{name}_{i}.png") for i in range(n_states)]
    widths, heights = zip(*(i.size for i in images))
    
    total_width = math.ceil(widths[0] * math.ceil(math.sqrt(n_states)))
    max_height = math.ceil(heights[0] * math.ceil(math.sqrt(n_states)))
    
    new_im = Image.new('RGB', (total_width, max_height))
    
    x_offset = 0
    y_offset = 0
    
    count = 0
    for im in images:
        count += 1
        
        new_im.paste(im, (x_offset, y_offset))
        x_offset += im.size[0]
        if (count == math.ceil(math.sqrt(n_states))):
            count = 0
            
            x_offset = 0
            y_offset += im.size[1]
            
    new_im.save(f"{plots}/{name}_Grid.png")

    
def gradientBoost(loss='squared_error', learning_rate=0.1, n_estimators=100, subsample=1.0, min_samples_split=2, min_samples_leaf=1, min_weight_fraction_leaf=0.0, max_depth=3, min_impurity_decrease=0.0, init=None, random_state=42, max_features=None, alpha=0.9, verbose=0, max_leaf_nodes=None, warm_start=False, validation_fraction=0.1, n_iter_no_change=None, tol=0.0001, ccp_alpha=0.0):
    local_parameters = locals()
    printToLog(f"# INFO - Running GBR with parameters [{local_parameters}]")
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=random_state)
    y_train = y_train.drop('REFCODE', axis=1)
    
    regressor = GradientBoostingRegressor(**local_parameters)
    if choices.__contains__("2"):
        regressor = MultiOutputRegressor(regressor)
    else:
        y_train = y_train.values.ravel()
         
    regressor.fit(X_train, y_train)
    
    y_pred = regressor.predict(X_test)
    
    mse = mean_squared_error(y_test.drop('REFCODE', axis=1), y_pred)
    r2 = r2_score(y_test.drop('REFCODE', axis=1), y_pred)
    
    prediction_data = os.path.join(prediction, f"_prediction_data_{random_state}.csv")
    if os.path.isfile(prediction_data):
        os.remove(prediction_data)
    printToLog("# INFO - Creating new sheet ["+ prediction_data + "]")

    with open(prediction_data, 'a') as file:
        headers = ["REFCODE"]
        for feature_predict in predict:
            headers.append(f"Actual {feature_predict}")
            headers.append(f"Predicted {feature_predict}")
            headers.append(f"Difference {feature_predict}")

        print(','.join(headers), file=file)
        
        for i in range(len(y_pred)):
            refcode = str(y_test.iloc[i]['REFCODE'])
            pred_df = pd.DataFrame(y_pred)

            data = [refcode]
            for feature_predict in predict:
                actual = float(y_test.iloc[i][feature_predict])
                predicted = pred_df[0][i]
                difference = abs(actual-predicted)            

                data.append(f"{actual:.2f}")
                data.append(f"{predicted:.2f}")
                data.append(f"{difference:.2f}")

            print(','.join(data), file=file)
    
        print(f"Parameters,\"{local_parameters}\"",file=file)
        print(f"Mean Squared Error,{mse:.2f}",file=file)
        print(f"R-squared Score,{r2:.2f}",file=file)
    

    perm_importance = permutation_importance(regressor,X_train,y_train,n_repeats=n_repeats,random_state=random_state)
    with open(mean_importance_summary, 'a') as file:
        data = [random_state]
        data.extend(list(perm_importance.importances_mean))

        print(','.join([str(x) for x in data]), file=file)
    with open(deviation_importance_summary, 'a') as file:
        data = [random_state]
        data.extend(list(perm_importance.importances_std))

        print(','.join([str(x) for x in data]), file=file)
    with open(prediction_summary, 'a') as file:
        print(f"{random_state},{mse:.2f},{r2:.2f},\"{local_parameters}\"", file=file)
       
    plotImportance(random_state)
    plotPrediction(random_state)
    joblib.dump(regressor, f"{models}/_model_{random_state}.pkl") 
    return (mse, r2)

def objective(trial, random_state=42):
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
    else:
        y_train = y_train.values.ravel()
        
    regressor.fit(X_train, y_train)
    y_pred = regressor.predict(X_test)
    return mean_squared_error(y_test.drop('REFCODE', axis=1), y_pred) 



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

log_keys = False
bond_cutoff = 2

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

feature_data = os.path.join(homeDirectory, "_feature_data.csv")
if choices.__contains__("1"):
    manual_inspect = os.path.join(homeDirectory, "manual_inspect")
    createDirectory(manual_inspect, "# WARN - No directory found for files to inspect manually.", False)
    
    if os.path.isfile(feature_data):
        os.remove(feature_data)
    printToLog("# INFO - Creating new sheet ["+ feature_data + "]")
    with open(feature_data, 'a') as file:
        print("REFCODE,W sigma,W sigma_11,W sigma_22,W sigma_33,P sigma,P sigma_11,P sigma_22,P sigma_33,W Distortion zeta,W Distortion delta,W Distortion sigma,W Distortion theta,C1-W-C2 angle,C1-W dist,C1-O dist,W-C1-O angle,C1 sigma,C1 sigma_11,C1 sigma_22,C1 sigma_33,C2-W dist,C2-O dist,W-C2-O angle,C2 sigma,C2 sigma_11,C2 sigma_22,C2 sigma_33,C3-W-C4 angle,C3-W dist,C3-O dist,W-C3-O angle,C3 sigma,C3 sigma_11,C3 sigma_22,C3 sigma_33,C4-W dist,C4-O dist,W-C4-O angle,C4 sigma,C4 sigma_11,C4 sigma_22,C4 sigma_33,C5-W-P angle,C5-W dist,C5-O dist,W-C5-O angle,C5 sigma,C5 sigma_11,C5 sigma_22,C5 sigma_33,P Distortion zeta,P Distortion delta,P Distortion sigma,P bonds,P-W dist,?1 mass,?1 bonds,W-P-?1 angle,P-?1 dist,?2 mass,?2 bonds,W-P-?2 angle,P-?2 dist,?3 mass,?3 bonds,W-P-?3 angle,P-?3 dist", file=file)
    
    output_files = os.path.join(os.path.join(homeDirectory, "MAIN"), "output_files")
    createDirectory(output_files, "# WARN - No directory found for output files.", True)
    directories = sorted([directory for directory in os.listdir(output_files) if os.path.isdir(os.path.join(output_files, directory))])
    
    if len(directories) == 0:
        printToLog("# WARN - No directories found in ["+ output_files + "]")
        quit()
    else:
        print(directories)
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
    
                                W_attached = getAttached(W)
                                symbols = list(map(lambda atom: atom.GetSymbol(), W_attached))                                
                                expected = ["P", "C", "C", "C", "C", "C"]
                                
                                printToLog(f"# INFO - Compound [{refcode}] Atom W [#{W_id}] has the following neighbours [{symbols}]")
                                if len(W_attached) == 6 and symbols == expected:   
                                    data[f"Site ID"] = f"{refcode}_{site_id}"
                                    data[f"W sigma"] = summary_line.split("(")[1].split(")")[0]
                                    data[f"W sigma_11"] = summary_line.split("[sigma_11 ")[1].split("]")[0]
                                    data[f"W sigma_22"] = summary_line.split("[sigma_22 ")[1].split("]")[0]
                                    data[f"W sigma_33"] = summary_line.split("[sigma_33 ")[1].split("]")[0]

                                    P = list(filter(lambda x: x.GetSymbol() == "P", W_attached))[0]
                                    
                                    for P_number, P_line in enumerate(summary_atoms, 0):
                                        if P_line.startswith(str(atoms[P.GetIdx()].strip().split()[-1])):
                                            data[f"P sigma"] = P_line.split("(")[1].split(")")[0]
                                            data[f"P sigma_11"] = P_line.split("[sigma_11 ")[1].split("]")[0]
                                            data[f"P sigma_22"] = P_line.split("[sigma_22 ")[1].split("]")[0]
                                            data[f"P sigma_33"] = P_line.split("[sigma_33 ")[1].split("]")[0]
                                            break

                                    W_site_coords = [parseCoord(conf.GetAtomPosition(W_id))]
                                    W_site_coords.extend(list(map(lambda atom: parseCoord(conf.GetAtomPosition(atom.GetIdx())), W_attached)))
                                    
                                    dist = oc.CalcDistortion(W_site_coords)
                                    data[f"W Distortion zeta"] = dist.zeta
                                    data[f"W Distortion delta"] = dist.delta
                                    data[f"W Distortion sigma"] = dist.sigma
                                    data[f"W Distortion theta"] = dist.theta

                                    carbons = sorted(W_attached[1:].copy(), key=lambda x: x.GetIdx() == getOpposite(conf, P, W).GetIdx())
                                    processed = []
                                    for C in carbons:                                    
                                        opposite = getOpposite(conf, C, W)
                                        if processed.__contains__(C.GetIdx()) or processed.__contains__(opposite.GetIdx()):
                                            continue

                                        processed.append(C.GetIdx())
                                        processed.append(opposite.GetIdx())

                                        # Angle between oppositeosite sides of W
                                        data[f"{C.GetSymbol()} {W.GetSymbol()} {opposite.GetSymbol()} [{C.GetIdx()}] [{W.GetIdx()}] [{opposite.GetIdx()}]"] = GetAngleDeg(conf,C.GetIdx(),W.GetIdx(),opposite.GetIdx())
    
                                        # Distance to W
                                        data[f"{C.GetSymbol()} {W.GetSymbol()} [{C.GetIdx()}] [{W.GetIdx()}]"] = GetBondLength(conf,C.GetIdx(),W.GetIdx())
                                        
                                        # Bonds to atom bonded to W
                                        for O in getAttached(C)[1:]:
                                            data[f"{C.GetSymbol()} {O.GetSymbol()} [{C.GetIdx()}] [{O.GetIdx()}]"] = GetBondLength(conf,C.GetIdx(),O.GetIdx())
                                            data[f"{W.GetSymbol()} {C.GetSymbol()} {O.GetSymbol()} [{W.GetIdx()}] [{C.GetIdx()}] [{O.GetIdx()}]"] = GetAngleDeg(conf,W.GetIdx(),C.GetIdx(),O.GetIdx())

                                        for C_number, C_line in enumerate(summary_atoms, 0):
                                            if C_line.startswith(str(atoms[C.GetIdx()].strip().split()[-1])):
                                                data[f"C [{C.GetIdx()}] sigma"] = C_line.split("(")[1].split(")")[0]
                                                data[f"C [{C.GetIdx()}] sigma_11"] = C_line.split("[sigma_11 ")[1].split("]")[0]
                                                data[f"C [{C.GetIdx()}] sigma_22"] = C_line.split("[sigma_22 ")[1].split("]")[0]
                                                data[f"C [{C.GetIdx()}] sigma_33"] = C_line.split("[sigma_33 ")[1].split("]")[0]
                                                break
                                        
                                        if not opposite.GetIdx() == P.GetIdx():
                                            data[f"{opposite.GetSymbol()} {W.GetSymbol()} [{opposite.GetIdx()}] [{W.GetIdx()}]"] = GetBondLength(conf,opposite.GetIdx(),W.GetIdx())
                                            for O in getAttached(opposite)[1:]:
                                                data[f"{opposite.GetSymbol()} {O.GetSymbol()} [{opposite.GetIdx()}] [{O.GetIdx()}]"] = GetBondLength(conf,opposite.GetIdx(),O.GetIdx())
                                                data[f"{W.GetSymbol()} {opposite.GetSymbol()} {O.GetSymbol()} [{W.GetIdx()}] [{opposite.GetIdx()}] [{O.GetIdx()}]"] = GetAngleDeg(conf,W.GetIdx(),opposite.GetIdx(),O.GetIdx())
                                            for C_number, C_line in enumerate(summary_atoms, 0):
                                                if C_line.startswith(str(atoms[opposite.GetIdx()].strip().split()[-1])):
                                                    data[f"C [{opposite.GetIdx()}] sigma"] = C_line.split("(")[1].split(")")[0]
                                                    data[f"C [{opposite.GetIdx()}] sigma_11"] = C_line.split("[sigma_11 ")[1].split("]")[0]
                                                    data[f"C [{opposite.GetIdx()}] sigma_22"] = C_line.split("[sigma_22 ")[1].split("]")[0]
                                                    data[f"C [{opposite.GetIdx()}] sigma_33"] = C_line.split("[sigma_33 ")[1].split("]")[0]
                                                    break
                                    
                                    P_attached = getAttached(P)
                                    symbols = list(map(lambda atom: atom.GetSymbol(), P_attached))                                
                                    printToLog(f"# INFO - Compound [{refcode}] Atom W [#{W_id}] attached P has the following neighbours [{symbols}]")

                                    if len(P_attached) <= bond_cutoff or len(P_attached) > 4:
                                        printToLog(f"# WARN - Compound [{refcode}] Atom P [{W_id}] has [{len(P_attached)}] P bonds")
                                        if not os.path.isfile(os.path.join(manual_inspect, refcode+'_cell.mol2')):
                                            shutil.copy(os.path.join(refcodeDirectory, refcode+'_cell.mol2'), manual_inspect)
                                        break
                                    P_site_dist = list(map(lambda atom: GetBondLength(conf,P.GetIdx(),atom.GetIdx()), P_attached))
                                    P_site_dist_mean = sum(dist for dist in P_site_dist) / len(P_site_dist)

                                    data[f"P Distortion zeta"] = sum(abs(dist - P_site_dist_mean) for dist in P_site_dist)
                                    data[f"P Distortion delta"] = sum(pow((dist - P_site_dist_mean) / P_site_dist_mean, 2) for dist in P_site_dist) / len(P_site_dist)
                                    P_site_angle = []
                                    for i in range(len(P_attached)):
                                        for j in range(i + 1, len(P_attached)):
                                            P_site_angle.append(GetAngleDeg(conf,P_attached[i].GetIdx(),P.GetIdx(),P_attached[j].GetIdx()))
                                    data[f"P Distortion sigma"] = sum(abs(109.5 - angle) for angle in P_site_angle)
          
                                    data[f"Bonds {P.GetSymbol()} {P.GetIdx()}"] = len(P_attached)
                                    data[f"{P.GetSymbol()} {W.GetSymbol()} [{P.GetIdx()}] [{W.GetIdx()}]"] = GetBondLength(conf,P.GetIdx(),W.GetIdx())

                                    # Bonds to P bonded to W
                                    for atom in P_attached:
                                        if not atom.GetIdx() == W.GetIdx():
                                            data[f"Mass {atom.GetSymbol()} {atom.GetIdx()}"] = atom.GetAtomicNum()
                                            data[f"Bonds {atom.GetSymbol()} {atom.GetIdx()}"] = int(len(getAttached(atom)))
                                            data[f"{W.GetSymbol()} {P.GetSymbol()} {atom.GetSymbol()} [{W.GetIdx()}] [{P.GetIdx()}] [{atom.GetIdx()}]"] = GetAngleDeg(conf,W.GetIdx(),P.GetIdx(),atom.GetIdx())
                                            data[f"{P.GetSymbol()} {atom.GetSymbol()} [{P.GetIdx()}] [{atom.GetIdx()}]"] = GetBondLength(conf,P.GetIdx(),atom.GetIdx())

                                    df = pd.read_csv(feature_data)
                                                
                                    values = list(data.values())
                                    while len(values) < len(df.columns):
                                        values.extend(["0","0","0","0"])
                                    df.loc[len(df)] = values


                                    if log_keys:
                                        keys = list(data.keys())
                                        while len(keys) < len(df.columns):
                                            keys.extend(["Placeholder Mass","Placeholder Bonds","Placeholder Angle","Placeholder Distance"])
                                        df.loc[len(df)] = keys
                                    df.to_csv(feature_data, index=False)
                                else:
                                    printToLog(f"# WARN - Compound [{refcode}] Atom W [{W_id}] invalid")
    
                    else:
                        printToLog(f"# WARN - Compound [{refcode}] Not complete")
            else:
                printToLog(f"# WARN - Compound [{refcode}] No summary file found")                            
if not os.path.isfile(feature_data):
    printToLog("# WARN - No .csv found.")
    quit()

predict = []

#Params
drop = ["P sigma","C1-W-C2 angle","C1-W dist","C1-O dist","W-C1-O angle","C1 sigma","C2-W dist","C2-O dist","W-C2-O angle","C2 sigma","C3-W-C4 angle","C3-W dist","C3-O dist","W-C3-O angle","C3 sigma","C4-W dist","C4-O dist","W-C4-O angle","C4 sigma","C5-W-P angle","C5-W dist","C5-O dist","W-C5-O angle","C5 sigma"]
#drop = ["P sigma","C1-W-C2 angle","C1-W dist","C1-O dist","W-C1-O angle","C1 sigma_11","C1 sigma_22","C1 sigma_33","C2-W dist","C2-O dist","W-C2-O angle","C2 sigma_11","C2 sigma_22","C2 sigma_33","C3-W-C4 angle","C3-W dist","C3-O dist","W-C3-O angle","C3 sigma_11","C3 sigma_22","C3 sigma_33","C4-W dist","C4-O dist","W-C4-O angle","C4 sigma_11","C4 sigma_22","C4 sigma_33","C5-W-P angle","C5-W dist","C5-O dist","W-C5-O angle","C5 sigma_11","C5 sigma_22","C5 sigma_33"]

#drop = ["P sigma","C1-W-C2 angle","C1-W dist","C1-O dist","W-C1-O angle","C1 sigma","C1 sigma_11","C1 sigma_22","C1 sigma_33","C2-W dist","C2-O dist","W-C2-O angle","C2 sigma","C2 sigma_11","C2 sigma_22","C2 sigma_33","C3-W-C4 angle","C3-W dist","C3-O dist","W-C3-O angle","C3 sigma","C3 sigma_11","C3 sigma_22","C3 sigma_33","C4-W dist","C4-O dist","W-C4-O angle","C4 sigma","C4 sigma_11","C4 sigma_22","C4 sigma_33","C5-W-P angle","C5-W dist","C5-O dist","W-C5-O angle","C5 sigma","C5 sigma_11","C5 sigma_22","C5 sigma_33"]
#drop = [] # ["?3 mass","?3 bonds","W-P-?3 angle","P-?3 dist"] #['W sigma','W sigma_11','W sigma_22','W sigma_33']

df = pd.read_csv(feature_data)
df = df.drop(drop, axis=1)

#Multiregressor target
if choices.__contains__("2"):
    df = df.drop(['W sigma'], axis=1)
    predict = ['W sigma_11','W sigma_22','W sigma_33']

#Single regressor target
elif choices.__contains__("3"):
    df = df.drop(['W sigma_11', 'W sigma_22', 'W sigma_33'], axis=1)
    predict = ['W sigma']

#Multiregressor target
#if choices.__contains__("2"):
#    df = df.drop(['P sigma'], axis=1)
#    predict = ['P sigma_11','P sigma_22','P sigma_33']

#Single regressor target
#elif choices.__contains__("3"):
#    df = df.drop(['P sigma_11', 'P sigma_22', 'P sigma_33'], axis=1)
#    predict = ['P sigma']

X = df.drop(['REFCODE'] + predict, axis=1)
y = df[['REFCODE'] + predict]

printToLog(X)
printToLog(y)
feature_names = list(X.columns.values)

if choices.__contains__("2") or choices.__contains__("3"):
    out_type = "single"
    if choices.__contains__("2"):
        out_type = "multi"
        
    time = str(datetime.datetime.now().strftime("[%Y-%m-%d_%H-%M-%S]"))
    prediction = os.path.join(homeDirectory, f"prediction_{out_type}_{time}")
    createDirectory(prediction, "# WARN - No directory found for prediction files.", False)

    plots = os.path.join(prediction, "plots")
    createDirectory(plots, "# WARN - No directory found for plots.", False)

    models = os.path.join(prediction, "models")
    createDirectory(models, "# WARN - No directory found for models.", False)

    plotCorrelation()
    X.to_csv(os.path.join(prediction,"data_X.csv"), index=False)
    y.to_csv(os.path.join(prediction,"data_y.csv"), index=False)

    prediction_summary = os.path.join(prediction, "prediction_summary.csv")
    if os.path.isfile(prediction_summary):
        os.remove(prediction_summary)
    with open(prediction_summary, 'a') as file:
        print(f"Random state,MSE,R-Squared,Parameters", file=file)


    headers = ["Random state"]
    headers.extend(feature_names)

    mean_importance_summary = os.path.join(prediction, "mean_importance_summary.csv")
    if os.path.isfile(mean_importance_summary):
        os.remove(mean_importance_summary)
    with open(mean_importance_summary, 'a') as file:
        print(','.join(headers), file=file)

    deviation_importance_summary = os.path.join(prediction, "deviation_importance_summary.csv")
    if os.path.isfile(deviation_importance_summary):
        os.remove(deviation_importance_summary)
    with open(deviation_importance_summary, 'a') as file:
        print(','.join(headers), file=file)
    
    overall = []
    for random_state in range(n_states):
        study = optuna.create_study(direction='minimize')
        study.optimize(lambda trial: objective(trial, random_state=random_state), n_trials=n_trials)
        printToLog(f"# INFO - Best parameters determined to be [{study.best_params}]")
        overall.append(gradientBoost(**study.best_params, random_state=random_state))        

    mse = sum(value[0] for value in overall) / len(overall)
    r2 = sum(value[1] for value in overall) / len(overall)
    
    printToLog(f"# INFO - Averaged Mean Squared Error: {mse:.2f}")
    printToLog(f"# INFO - Averaged R-squared Score: {r2:.2f}")
        
    with open(prediction_summary, 'a') as file:
        print(f"Average,{mse:.2f},{r2:.2f}", file=file)

    df = pd.read_csv(mean_importance_summary, index_col=["Random state"])
    df = pd.concat([df, df.mean().to_frame('Average').T])
    df.to_csv(mean_importance_summary,index_label=["Random state"])

    df = pd.read_csv(deviation_importance_summary, index_col=["Random state"])
    df = pd.concat([df, df.mean().to_frame('Average').T])
    df.to_csv(deviation_importance_summary,index_label=["Random state"])

    plotPredictionSummary()
    plotImportanceSummary()
    createGrid("_prediction_data")
    createGrid("_importance_data")
