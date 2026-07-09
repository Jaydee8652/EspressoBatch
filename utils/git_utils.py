# Utility functions for reading and writing to a global csv stored on github

#Imports
from github import Auth
from github import Github
from github import GithubIntegration
from github import Github
import __main__ as main
import os
import subprocess
import sys
import re 
import pandas as pd
import io
import csv
import datetime
import time
from utils.generic_utils import *

# Params - can be changed
location = "Rorqual"

REPO = 'REPO_NAME' #eg 'ExampleRepo' for https://github.com/Jaydee8652/ExampleRepo 
sheetPath = 'sheet.csv'
flagPath = 'sheet_flag.txt'
auth = Auth.Token("github_pat_0000000000000000000000000000000000000000000000000000000000000000000000000000000000") #Must have permissions on repo 

# Git Authentication
g = Github(auth=auth)
g.get_user().login

homeDirectory = os.getcwd()#Directory where we are
localSheet = os.path.join(homeDirectory, sheetPath)
localFlag = os.path.join(homeDirectory, flagPath)

repo = g.get_user().get_repo(REPO)
all_files = []
contents = repo.get_contents("")
while contents:
    file_content = contents.pop(0)
    if file_content.type == "dir":
        contents.extend(repo.get_contents(file_content.path))
    else:
        file = file_content
        all_files.append(str(file).replace('ContentFile(path="','').replace('")',''))

# Set the flag on github
def setFlag(log, boolean):
    source = log.split(".")[0]

    if flagPath in all_files:
        flag = repo.get_contents(flagPath)
        if os.path.isfile(localFlag):
            os.remove(localFlag)        
        with open(localFlag, "w") as file:
            print(boolean, file=file)
        with open(localFlag, "r") as file:
            repo.update_file(flag.path, "AC at ["+str(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))+"] by "+str(source)+" ["+str(boolean)+"]", file.read(), flag.sha)
        os.remove(localFlag)        

#Download csv from github
def downloadCSV(log):
    gitContent = repo.get_contents(sheetPath).decoded_content.decode()

    if os.path.isfile(localSheet):
        printToLog(log, "# INFO - Removing existing local file ["+ localSheet + "]")
        os.remove(localSheet)# Clear current local copy
    printToLog(log, "# INFO - Downloading ["+sheetPath+"] at ["+sheetPath+"] from [REPO - "+REPO+"]")
    with open(localSheet, 'a') as file:
        file.write(gitContent)# Save data to local copy
    return localSheet

#Upload csv to github
def uploadCSV(log):
    git = repo.get_contents(sheetPath)

    with open(localSheet, 'r') as file:
        printToLog(log, "# INFO - Attempting to update ["+sheetPath+"] at ["+sheetPath+"] in [REPO - "+REPO+"]")
        source = log.split(".")[0]

        localContent = file.read()
        repo.update_file(git.path, "AC at ["+str(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))+"] by "+str(source)+" [sheet.csv]", localContent, git.sha)#Update global .csv
        setFlag(log, "True")
        printToLog(log, "# INFO - Updated ["+sheetPath+"] at ["+sheetPath+"] in [REPO - "+REPO+"]")
        if os.path.isfile(localSheet):
            os.remove(localSheet)# Clear current local copy

#Get current location, defined in gitutils
def getLocation():
    return location

# Reference the flag on github, ensures the global .csv is not altered by two scripts at once
def verify(log):
    if sheetPath in all_files and flagPath in all_files:
        printToLog(log, "# INFO - Requesting master .csv")
        flag = repo.get_contents(flagPath)
        flagContent = repo.get_contents(flagPath).decoded_content.decode()
        if(flagContent.strip() == "True"):
            printToLog(log, "# INFO - Master .csv availabile")
            setFlag(log, "False")
            return True
        else:
            printToLog(log, "# INFO - Master .csv currently in use. Waiting for availability")
            time.sleep(30)
            return verify(log)

#Create sheet
def initSheet(log):
    if os.path.isfile(localSheet):
        printToLog(log,"# INFO - Removing existing local file ["+ localSheet + "]")
        os.remove(localSheet)# Clear current local copy
            
        printToLog(log,"# INFO - Creating new sheet ["+ localSheet + "]")
        with open(localSheet, 'a') as file:
            file.write("[REFCODE]")
        df = pd.read_csv(localSheet)
        df = pd.concat([df, pd.DataFrame({"[REFCODE]": ["init"]})], ignore_index=True)    
    
        df["[BATCH_location]"] = ["Abyss"]    
        df["[BATCH_done]"] = ["True"]
        df["[BATCH_time]"] = ["2026-06-26 19:10:06"]
    
        df["[PWSCF_done]"] = ["True"]
        df["[PWSCF_time]"] = ["2026-06-26 19:10:06"]
    
        df["[GIPAW_done]"] = ["True"]
        df["[GIPAW_time]"] = ["2026-06-26 19:10:06"]
    
        df["[PWSCF_version]"] = ["v.7.3.1"]
        df["[PWSCF_numberMPI]"] = ["1"]
        df["[PWSCF_numberThreads]"] = ["1"]
        df["[PWSCF_RG]"] = ["1"]
        df["[PWSCF_estimatedRAM]"] = ["1.1"]
        df["[PWSCF_scfCycles]"] = ["1"]
        df["[PWSCF_bfgsSteps]"] = ["1"]
        df["[PWSCF_finalEnergy]"] = ["-1"]
        
        df["[GIPAW_version]"] = ["v.7.3.1"]
        df["[GIPAW_numberMPI]"] = ["1"]
        df["[GIPAW_numberThreads]"] = ["1"]
        df["[GIPAW_RG]"] = ["1"]
        df["[GIPAW_mscPPM]"] = ["1.1"]
        df["[GIPAW_msCorrection]"] = ["[['0.6667', '0.0000', '0.0000'], ['0.0000', '0.6667', '0.0000'], ['0.0000', '0.0000', '0.6667']]"]
        df.to_csv(localSheet, index=False)

#localLog = "utils.log"
#initSheet(localLog)
#uploadCSV(localLog)