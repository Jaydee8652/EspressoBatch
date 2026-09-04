# Trains a GradientBoostingRegressor model, saved in the "models" directory. 
# Called by create_model.py through slurm, do not run manually
import re 
import os
import sys
import time
import datetime
import pandas as pd
import numpy as np
import shutil
import math

#Jank thing to fix the path. Very annoying artefact of running python scripts by absolute path.
sys.path[0] = sys.path[0][:-6] + sys.path[0][-6:].replace("/utils", "")

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from matplotlib import colormaps
from PIL import Image

from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.svm import SVR
from sklearn.inspection import permutation_importance
from sklearn.multioutput import MultiOutputRegressor

from utils.generic_utils import printToLog as pl, createDirectory as cd, cellVolume, featureExtractor
from utils.params import *

import xgboost as xgb
import optuna
import joblib

n_trials = 100 # 100 - Number of parameter combinations trialled per seed
n_states = 100 # 100 - Number of seeds
n_repeats = 150 # 150 - Permutations for feature importance testing
datapoint_cap = 10000 # Cap number of datapoints to test on to this number or lower
consider = ["C5 sigma_33", "W Distortion zeta","W Distortion delta","W Distortion sigma","W Distortion theta","P Distortion zeta","P Distortion delta","P Distortion sigma","P bonds","P-W dist","?1 atomic number","?1 bonds","W-P-?1 angle","P-?1 dist","?2 atomic number","?2 bonds","W-P-?2 angle","P-?2 dist","?3 atomic number","?3 bonds","W-P-?3 angle","P-?3 dist"]

cmap = plt.cm.viridis

multiregressor = False
predictP = True # If false will instead predict W

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)

def plotPrediction(random_state):
    title = f"_prediction_data_{random_state}"
    printToLog(f"# INFO - Plotting [{title}]")

    plt.clf()
    prediction_data = os.path.join(directory, f"{title}.csv")
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
    plt.axline((min(xmin, ymin),min(xmin, ymin)), slope=1, linewidth=1, color='red')
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
    plt.title(f"{title} - Mean {np.mean(y)}")
    
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
    
    total_width = math.ceil(max([widths[i] for i in range(len(widths))]) * math.ceil(math.sqrt(n_states)))
    max_height = math.ceil(min([heights[i] for i in range(len(heights))]) * math.ceil(math.sqrt(n_states)))
    
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


def xgBoost(random_state=42, n_estimators=20):
    local_parameters = locals()
    printToLog(f"# INFO - Running XGB with parameters [{local_parameters}]")
    xgb.set_config(verbosity=2)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=random_state)
    y_train = y_train.drop('REFCODE', axis=1)

    regressor = xgb.XGBRegressor(**local_parameters)
    regressor.fit(X_train, y_train, eval_set=[(X_test, y_test.drop('REFCODE', axis=1))])

    y_pred = regressor.predict(X_test)
    
    mse = mean_squared_error(y_test.drop('REFCODE', axis=1), y_pred)
    r2 = r2_score(y_test.drop('REFCODE', axis=1), y_pred)
    
    prediction_data = os.path.join(directory, f"_prediction_data_{random_state}.csv")
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
        pred_df = pd.DataFrame(y_pred)

        for i in range(len(y_pred)):
            refcode = str(y_test.iloc[i]['REFCODE'])

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

    print(regressor.feature_importances_)

    with open(mean_importance_summary, 'a') as file:
        data = [random_state]
        data.extend(list(regressor.feature_importances_))

        print(','.join([str(x) for x in data]), file=file)
    with open(deviation_importance_summary, 'a') as file:
        data = [random_state]
        data.extend(len(list(regressor.feature_importances_)) * [0])

        print(','.join([str(x) for x in data]), file=file)
    with open(prediction_summary, 'a') as file:
        print(f"{random_state},{mse:.2f},{r2:.2f},\"{local_parameters}\"", file=file)
       
    plotImportance(random_state)
    plotPrediction(random_state)
    regressor.save_model(f"{models}/_model_{random_state}.json")
    return (mse, r2)

def xgb_objective(trial, random_state=42):
    param_space = {
        'n_estimators': trial.suggest_int('n_estimators', 50, 500, step=25),
        'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.2),
        'num_parallel_tree': trial.suggest_int('num_parallel_tree', 1, 6),
    }
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=random_state)
    y_train = y_train.drop('REFCODE', axis=1)

    regressor = xgb.XGBRegressor(**param_space, random_state=random_state)
    regressor.fit(X_train, y_train, eval_set=[(X_test, y_test.drop('REFCODE', axis=1))])
    y_pred = regressor.predict(X_test)

    return mean_squared_error(y_test.drop('REFCODE', axis=1), y_pred) 

def gradientBoost(loss='squared_error', learning_rate=0.1, n_estimators=100, subsample=1.0, min_samples_split=2, min_samples_leaf=1, min_weight_fraction_leaf=0.0, max_depth=3, min_impurity_decrease=0.0, init=None, random_state=42, max_features=None, alpha=0.9, verbose=0, max_leaf_nodes=None, warm_start=False, validation_fraction=0.1, n_iter_no_change=None, tol=0.0001, ccp_alpha=0.0):
    local_parameters = locals()
    printToLog(f"# INFO - Running GBR with parameters [{local_parameters}]")
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=random_state)
    y_train = y_train.drop('REFCODE', axis=1)
    
    regressor = GradientBoostingRegressor(**local_parameters)
    if multiregressor:
        regressor = MultiOutputRegressor(regressor)
    else:
        y_train = y_train.values.ravel()
         
    regressor.fit(X_train, y_train)
    
    y_pred = regressor.predict(X_test)

    mse = mean_squared_error(y_test.drop('REFCODE', axis=1), y_pred)
    r2 = r2_score(y_test.drop('REFCODE', axis=1), y_pred)
    
    prediction_data = os.path.join(directory, f"_prediction_data_{random_state}.csv")
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

def gb_objective(trial, random_state=42):
    param_space = {
        'n_estimators': trial.suggest_int('n_estimators', 50, 500, step=25),
        'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.2),
        'max_depth': trial.suggest_int('max_depth', 1, 6),
    }
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=random_state)
    y_train = y_train.drop('REFCODE', axis=1)

    regressor = GradientBoostingRegressor(**param_space, random_state=random_state)
    if multiregressor:
        regressor = MultiOutputRegressor(regressor)
    else:
        y_train = y_train.values.ravel()
        
    regressor.fit(X_train, y_train)
    y_pred = regressor.predict(X_test)
    return mean_squared_error(y_test.drop('REFCODE', axis=1), y_pred) 

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
directory = os.getcwd() #Directory where we are
homeDirectory = os.path.split(os.path.split(directory)[0])[0] 
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ directory + "]")    

createDirectory(directory, "# WARN - No directory found", True)
feature_data = os.path.join(homeDirectory, "_training_feature_data.csv")                     
if not os.path.isfile(feature_data):
    printToLog("# WARN - No .csv found.")
    quit()

predict = []
df = pd.read_csv(feature_data)
if predictP:
    df = df.filter(consider + ["REFCODE","P sigma","P sigma_11","P sigma_22","P sigma_33"])
    
    if multiregressor:
        df = df.drop(['P sigma'], axis=1)
        predict = ['P sigma_11','P sigma_22','P sigma_33']
    else:
        df = df.drop(['P sigma_11', 'P sigma_22', 'P sigma_33'], axis=1)
        predict = ['P sigma']
else:
    df = df.filter(consider + ["REFCODE","W sigma","W sigma_11","W sigma_22","W sigma_33"])
    
    if multiregressor:
        df = df.drop(['W sigma'], axis=1)
        predict = ['W sigma_11','W sigma_22','W sigma_33']
    else:
        df = df.drop(['W sigma_11', 'W sigma_22', 'W sigma_33'], axis=1)
        predict = ['W sigma']    


df = df.head(datapoint_cap)
X = df.drop(['REFCODE'] + predict, axis=1)
y = df[['REFCODE'] + predict]

printToLog(X)
printToLog(y)
feature_names = list(X.columns.values)

plots = os.path.join(directory, "plots")
createDirectory(plots, "# WARN - No directory found for plots.", False)

models = os.path.join(directory, "models")
createDirectory(models, "# WARN - No directory found for models.", False)

plotCorrelation()
X.to_csv(os.path.join(directory,"data_X.csv"), index=False)
y.to_csv(os.path.join(directory,"data_y.csv"), index=False)

prediction_summary = os.path.join(directory, "prediction_summary.csv")
if os.path.isfile(prediction_summary):
    os.remove(prediction_summary)
with open(prediction_summary, 'a') as file:
    print(f"Random state,MSE,R-Squared,Parameters", file=file)

headers = ["Random state"]
headers.extend(feature_names)

mean_importance_summary = os.path.join(directory, "mean_importance_summary.csv")
if os.path.isfile(mean_importance_summary):
    os.remove(mean_importance_summary)
with open(mean_importance_summary, 'a') as file:
    print(','.join(headers), file=file)

deviation_importance_summary = os.path.join(directory, "deviation_importance_summary.csv")
if os.path.isfile(deviation_importance_summary):
    os.remove(deviation_importance_summary)
with open(deviation_importance_summary, 'a') as file:
    print(','.join(headers), file=file)

overall = []
for random_state in range(n_states):
    study = optuna.create_study(direction='minimize')
    study.optimize(lambda trial: gb_objective(trial, random_state=random_state), n_trials=n_trials)
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