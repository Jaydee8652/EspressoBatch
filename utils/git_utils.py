# Utility functions for reading and writing to a global csv stored on github

#Imports
import __main__ as main
import os
import subprocess
import re 
import pandas as pd
import io
import csv
import datetime
import time
import sys
from utils.generic_utils import *

# Params - can be changed
REPO = "REPO_NAME" #https://github.com/Jaydee8652/REPO_NAME
TOKEN = "github_pat_0000000000000000000000000000000000000000000000000000000000000000000000000000000000" #Must have permissions on repo

sheetPath = 'sheet.csv'
flagPath = 'sheet_flag.txt'

#Main
homeDirectory = os.getcwd()#Directory where we are
localSheet = os.path.join(homeDirectory, sheetPath)
localFlag = os.path.join(homeDirectory, flagPath)

GIT_ACTIVE = False

# Git Authentication
if not TOKEN == "github_pat_0000000000000000000000000000000000000000000000000000000000000000000000000000000000" and not REPO == "REPO_NAME":
    #Imports
    from github import Auth
    from github import Github
    
    auth = Auth.Token(TOKEN)
    g = Github(auth=auth)
    g.get_user().login
    
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
    GIT_ACTIVE = True

local_log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
utils = os.path.join(homeDirectory, "utils")
location = os.path.join(utils, "location.txt")
if not os.path.exists(location): 
    with open(location, "a") as file:
        printToLog(local_log," --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - No location data found. Attempting to retreive.")    
        try:
            out = subprocess.check_output(['hostname'],shell=True)
            out = out.decode("utf-8").strip()
    
            print(out,file=file)
            printToLog(local_log,"# INFO - Location determined to be ["+str(out)+"], saved to ["+str(location)+"]")
            printToLog(local_log,"# INFO - Override manually by changing the contents of ['location.txt']")
        except subprocess.CalledProcessError as e:
            printToLog(local_log,"# INFO - Error retreiving llocation data.")
            printToLog(local_log,str(e))

def getLocation():
    with open(location, "r") as file:
        return file.read().strip()

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
    if GIT_ACTIVE:
        gitContent = repo.get_contents(sheetPath).decoded_content.decode()
    
        if os.path.isfile(localSheet):
            printToLog(log, "# INFO - Removing existing local file ["+ localSheet + "]")
            os.remove(localSheet)# Clear current local copy
        printToLog(log, "# INFO - Downloading ["+sheetPath+"] at ["+sheetPath+"] from [REPO - "+REPO+"] - DO NOT CANCEL")
        with open(localSheet, 'a') as file:
            file.write(gitContent)# Save data to local copy
        return localSheet
    else:
        printToLog(log, "# INFO - Git integration inactive. Retreiving local .csv ["+sheetPath+"]")
        return localSheet



#Upload csv to github
def uploadCSV(log):
    if GIT_ACTIVE:
        git = repo.get_contents(sheetPath)
    
        with open(localSheet, 'r') as file:
            printToLog(log, "# INFO - Attempting to update ["+sheetPath+"] at ["+sheetPath+"] in [REPO - "+REPO+"] - DO NOT CANCEL")
            source = log.split(".")[0]
    
            localContent = file.read()
            repo.update_file(git.path, "AC at ["+str(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))+"] by "+str(source)+" [sheet.csv]", localContent, git.sha)#Update global .csv
            setFlag(log, "True")
            printToLog(log, "# INFO - Updated ["+sheetPath+"] at ["+sheetPath+"] in [REPO - "+REPO+"]")
            if os.path.isfile(localSheet):
                os.remove(localSheet)# Clear current local copy
    else:
        printToLog(log, "# INFO - Git integration inactive. Updated local .csv ["+sheetPath+"]")


    
# Reference the flag on github, ensures the global .csv is not altered by two scripts at once
def verify(log):
    if GIT_ACTIVE:
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
    else:
        printToLog(log, "# INFO - Git integration inactive. Attempting to retrieve local .csv")
        if not os.path.isfile(localSheet):
            printToLog(log, "# INFO - No local .csv found. Attempting to initialise")
            initSheet(log)
        return True



# Appends the refcode of all local input directories to a local .csv to be referenced by other scripts
def appendCSV(log):
    processedCount = 0    
    printToLog(log,"# INFO - Attempting to append sheet at ["+str(localSheet)+"]")

    if not os.path.isfile(localSheet):
        printToLog(log, "# WARN - No local .csv found.")
        quit()
        
    inputPath = os.path.join(homeDirectory, "Input_Files")
    createDirectory(log, inputPath, "# WARN - No directory found for input files.", True)
    directories = [directory for directory in os.listdir(inputPath) if os.path.isdir(os.path.join(inputPath, directory)) and not directory.startswith(".") and not os.path.isfile(os.path.join(os.path.join(inputPath, directory), "INCOMPLETE.txt"))]
    
    printToLog(log,"# INFO - The following input directories are available ["+str(directories)+"]")
    
    df = pd.read_csv(localSheet)
    for refcode in directories:
        if refcode in df['[REFCODE]'].values:
            printToLog(log,"# INFO - Compound ["+ refcode +"] Already present in sheet")
        else:           
            printToLog(log,"# INFO - Compound ["+ refcode +"] Appending to sheet")
            df = pd.concat([df, pd.DataFrame({"[REFCODE]": [refcode]})], ignore_index=True)
            processedCount += 1

    df.to_csv(localSheet, index=False)
    printToLog(log,"# INFO - Successfully appended ["+str(processedCount)+"] compounds to sheet at ["+str(localSheet)+"]")



# Extracts data summary file and updates a local .csv
def updateCSV(log):
    processedCount = 0
    printToLog(log,"# INFO - Attempting to update sheet at ["+str(localSheet)+"]")
    if not os.path.isfile(localSheet):
        printToLog(log, "# WARN - No local .csv found.")
        quit()
    
    summariesPath = os.path.join(homeDirectory, "Summary_Files")
    createDirectory(log,summariesPath, "# INFO - No directory found for summary files, creating at ["+str(summariesPath)+"]", False)
    summaryFiles = [file for file in os.listdir(summariesPath) if file.endswith('_summary.txt') and os.path.isfile(os.path.join(summariesPath, file))]#Get .UPFs from directory
    printToLog(log,"# INFO - The following summaries are available ["+str(summaryFiles)+"]")

    df = pd.read_csv(localSheet)  
    df.set_index('[REFCODE]', inplace = True)
    df = df.astype(str)
    
    for summary in summaryFiles:
        refcode = os.path.splitext(summary)[0].replace("_summary", "")
        printToLog(log,"# INFO - Compound ["+ refcode +"] Processing output data")

        with open(os.path.join(summariesPath, summary), "r") as file:            
            read = file.read()
            lines = read.splitlines()
            
            for line in lines:
                if not len(line) == 0 and not line.startswith("#") and not line.startswith("_"):
                    value = line[line.find("=")+1:].strip()
                    name = line[:line.find("=")-1].strip()
                    writeCSV(df, refcode, "["+str(name)+"]", str(value))
            processedCount += 1

    df = df.replace("nan", "")
    df.to_csv(localSheet)#Update local csv
    printToLog(log,"# INFO - Successfully updated data in sheet at ["+str(localSheet)+"] for ["+str(processedCount)+"] compounds")



# References and updates a local .csv to submit requests to slurm, only running calculations not already flagged as batched
# Batches 'batchCount' every run to avoid requesting too many resources at once 
def batchCalculations(log, batchCount):
    printToLog(log,"# INFO - Attempting to batch ["+str(batchCount)+"] calculations")
    if not os.path.isfile(localSheet):
        printToLog(log, "# WARN - No local .csv found.")
        quit()
    
    #Make sure there is a directory to process
    inputPath = os.path.join(homeDirectory, "Input_Files")
    createDirectory(log,inputPath, "# WARN - No directory found for input files.", True)
    directories = [directory for directory in os.listdir(inputPath) if os.path.isdir(os.path.join(inputPath, directory)) and not directory.startswith(".") and not os.path.isfile(os.path.join(os.path.join(inputPath, directory), "INCOMPLETE.txt"))]

    numberOfDirectories = len(directories) # determine number of directories
    if numberOfDirectories == 0:
        printToLog(log,"# WARN - No directories found in ["+ inputPath + "]")
        quit()
    else:
        printToLog(log,"# INFO - [" + str(numberOfDirectories) + "] directories found at ["+ inputPath + "]")
        
    processedCount = 0

    df = pd.read_csv(localSheet)
    df.set_index('[REFCODE]', inplace = True)

    for refcode in directories:
        if processedCount < batchCount:
            printToLog(log,"# INFO - Processing compound with refcode ["+ refcode +"]")
            if not df.at[refcode, "[BATCH_done]"] == "True":
                printToLog(log,"# INFO - Compound with refcode ["+ refcode +"] not previously run, attempting to batch")
                refcodeDirectory = os.path.join(inputPath, refcode)
                
                QE_SUB = os.path.join(refcodeDirectory, "QE_SUB")
                batch_path = os.path.join(refcodeDirectory, refcode+"_batch.txt")

                batchCommand = f"module load {getModules()}; cd {refcodeDirectory}; sbatch QE_SUB"
                if os.path.exists(QE_SUB):
                    try:
                        now = str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                        subprocess.call(batchCommand,shell=True)
                        
                        with open(batch_path, "a") as batch:
                            writeCSV(df, refcode, "[BATCH_done]", True)
                            writeCSV(df, refcode, "[BATCH_time]", now)
                            writeCSV(df, refcode, "[BATCH_location]", getLocation())
                            
                            print("\n# -Batch data\n", file=batch)
                            print("BATCH_done = "+str(True), file=batch)
                            print("BATCH_time = "+str(now), file=batch)
                            print("BATCH_location = "+str(getLocation()), file=batch)
                            
                            printToLog(log,"# INFO - Successfully batched calculation for compound ["+refcode+"] at ["+str(now)+"] on ["+str(getLocation())+"]")
                            processedCount += 1
                    except subprocess.CalledProcessError as e:
                        printToLog(log,"# WARN - Error batching calculation for compound with refcode ["+refcode+"]")
                        printToLog(str(e))
                else:
                    printToLog(log,"# WARN - QE_SUB not present for compound with refcode ["+refcode+"]")
            else:
                printToLog(log,"# INFO - Compound with refcode ["+refcode+"] has been previously batched at ["+str(df.at[refcode, "[BATCH_time]"])+"] on ["+str(df.at[refcode, "[BATCH_location]"])+"]")  
    df.to_csv(localSheet)#Update local csv
    printToLog(log,"# INFO - ["+str(processedCount)+"] Calculations successfully batched.")
    if processedCount < batchCount:
        printToLog(log,"# INFO - No more calculations to batch!")
    getQueue(log)


def getQueue(log):
    printToLog(log,"# INFO - Attempting to retrieve current slurm queue.")
    length = 0
    try:
        out = subprocess.check_output(['squeue --me'],shell=True)
        out = out.decode("utf-8")
    
        lines = out.splitlines()
        for line in lines:
            if "_SUB" in line:
                length += 1
            printToLog(log, line)
        printToLog(log,"# INFO - Slurm queue contains ["+str(length)+"] batched calculations.")
        return length
    except subprocess.CalledProcessError as e:
        printToLog(log,"# INFO - Error retreiving slurm queue.")
        printToLog(log,str(e))

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
    df["[PWSCF_done]"] = ["True"]
    df["[GIPAW_done]"] = ["True"]

    df["[BATCH_time]"] = ["2026-06-26 19:10:06"]
    df["[PWSCF_time]"] = ["2026-06-26 19:10:06"]
    df["[GIPAW_time]"] = ["2026-06-26 19:10:06"]

    df["[CIF_symmetryEquivalents]"] = ["2"]
    df["[CIF_symmetryFactor]"] = ["2"]
    df["[CIF_inversionCentre]"] = ["True"]

    df["[PWSCF_ecutwfc]"] = ["55.0"]
    df["[PWSCF_ecutrho]"] = ["8"]
    df["[PWSCF_conv_thr]"] = ["1.D-6"]
    df["[PWSCF_version]"] = ["v.7.3.1"]
    df["[PWSCF_numberMPI]"] = ["1"]
    df["[PWSCF_numberThreads]"] = ["1"]
    df["[PWSCF_RG]"] = ["1"]
    df["[PWSCF_estimatedRAM]"] = ["1.1"]
    df["[PWSCF_scfCycles]"] = ["1"]
    df["[PWSCF_bfgsSteps]"] = ["1"]
    df["[PWSCF_finalEnergy]"] = ["-1"]
    
    df["[GIPAW_q_gipaw]"] = ["0.01"]
    df["[GIPAW_version]"] = ["v.7.3.1"]
    df["[GIPAW_numberMPI]"] = ["1"]
    df["[GIPAW_numberThreads]"] = ["1"]
    df["[GIPAW_RG]"] = ["1"]
    df["[GIPAW_mscPPM]"] = ["1.1"]
    df["[GIPAW_msCorrection]"] = ["[['0.6667', '0.0000', '0.0000'], ['0.0000', '0.6667', '0.0000'], ['0.0000', '0.0000', '0.6667']]"]
    df.to_csv(localSheet, index=False)



