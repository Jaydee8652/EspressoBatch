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


#Jank thing to fix the path. Very annoying artefact of running python scripts by absolute path.
sys.path[0] = sys.path[0][:-6] + sys.path[0][-6:].replace("/utils", "")

from utils.generic_utils import *
from utils.params import *

#Main
homeDirectory = os.getcwd() # Directory where we are

qe_params = os.path.join(homeDirectory, "_qe_params.csv")
localSheet = os.path.join(homeDirectory, param_sheetPath)
localFlag = os.path.join(homeDirectory, param_flagPath)

GIT_ACTIVE = False

# Git Authentication - if the params for github have been changed, try to log in with them 
if not param_token == "github_pat_0000000000000000000000000000000000000000000000000000000000000000000000000000000000" and not param_repo == "REPO_NAME":
    #Imports
    from github import Auth
    from github import Github
    
    auth = Auth.Token(param_token)
    g = Github(auth=auth)
    g.get_user().login
    
    repo = g.get_user().get_repo(param_repo)
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

# Set the flag on github = boolean True/False with True being "is available"
def setFlag(log, boolean):
    source = log.split(".")[0]

    if param_flagPath in all_files:
        flag = repo.get_contents(param_flagPath)
        if os.path.isfile(localFlag):
            os.remove(localFlag)        
        with open(localFlag, "w") as file:
            print(boolean, file=file)
        with open(localFlag, "r") as file:
            repo.update_file(flag.path, "AC at ["+str(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))+"] by "+str(source)+" on "+str(param_location)+" ["+str(boolean)+"]", file.read(), flag.sha)
        os.remove(localFlag)        

# Reference the flag on github, ensures the global .csv is not altered by two scripts at once
def verify(log):
    if GIT_ACTIVE:
        if param_sheetPath in all_files and param_flagPath in all_files:
            printToLog(log, "# INFO - Requesting master .csv")
            flag = repo.get_contents(param_flagPath)
            flagContent = repo.get_contents(param_flagPath).decoded_content.decode()
            if(flagContent.strip() == "True"):
                printToLog(log, "# INFO - Master .csv availabile")
                setFlag(log, "False")
                return True
            else:
                printToLog(log, "# INFO - Master .csv currently in use. Waiting for availability")
                time.sleep(30) # Wait 30 seconds before trying again
                return verify(log)
    else:
        printToLog(log, "# INFO - Git integration inactive. Attempting to retrieve local .csv")
        if not os.path.isfile(localSheet):
            printToLog(log, "# INFO - No local .csv found. Attempting to initialise")
            initSheet(log)
        return True

#Download csv from github
def downloadCSV(log):
    if GIT_ACTIVE:
        gitContent = repo.get_contents(param_sheetPath).decoded_content.decode()
    
        if os.path.isfile(localSheet):
            printToLog(log, "# INFO - Removing existing local file ["+ localSheet + "]")
            os.remove(localSheet) # Clear current local copy
        printToLog(log, "# INFO - Downloading ["+param_sheetPath+"] at ["+param_sheetPath+"] from [REPO - "+param_repo+"] - DO NOT CANCEL")
        with open(localSheet, 'a') as file:
            file.write(gitContent) # Save data to local copy
        return localSheet
    else:
        printToLog(log, "# INFO - Git integration inactive. Retreiving local .csv ["+param_sheetPath+"]")
        return localSheet

#Upload csv to github
def uploadCSV(log):
    if GIT_ACTIVE:
        git = repo.get_contents(param_sheetPath)
    
        with open(localSheet, 'r') as file:
            printToLog(log, "# INFO - Attempting to update ["+param_sheetPath+"] at ["+param_sheetPath+"] in [REPO - "+param_repo+"] - DO NOT CANCEL")
            source = log.split(".")[0]
    
            localContent = file.read()
            repo.update_file(git.path, "AC at ["+str(str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))+"] by "+str(source)+" on "+str(param_location)+" [sheet.csv]", localContent, git.sha) # Update global .csv
            setFlag(log, "True")
            printToLog(log, "# INFO - Updated ["+param_sheetPath+"] at ["+param_sheetPath+"] in [REPO - "+param_repo+"]")
            if os.path.isfile(localSheet):
                os.remove(localSheet) # Clear current local copy
    else:
        printToLog(log, "# INFO - Git integration inactive. Updated local .csv ["+param_sheetPath+"]")

# Appends the refcode of all validated cif files to a local .csv to be referenced by other scripts 
# if github is enables the local .csv is a copy of the global .csv that will replace it
def appendCSV(log):
    processedCount = 0    
    printToLog(log,"# INFO - Attempting to append sheet at ["+str(localSheet)+"]")

    if not os.path.isfile(localSheet):
        printToLog(log, "# WARN - No local .csv found.")
        quit()

    if not os.path.isfile(qe_params):
        printToLog(log,"# WARN - No qe_params.csv found.")
        quit()
    
    cifs_path = os.path.join(homeDirectory, "cifs")
    createDirectory(log,cifs_path, "# WARN - No directory found for .cifs to process.", False)
    
    validated = os.path.join(cifs_path, "validated")
    createDirectory(log,validated, "# WARN - No directory found for .cifs to process. Place .cif files in or replace the newly created directory at ["+validated+"]", True)

    cifs = sorted([os.path.splitext(file)[0].replace(".cif", "") for file in sorted(os.listdir(validated)) if file.endswith('.cif') and os.path.isfile(os.path.join(validated, file))])

    printToLog(log,"# INFO - The following input directories are available ["+str(cifs)+"]")
    
    df = pd.read_csv(localSheet)

    qe = pd.read_csv(qe_params, encoding="utf-8-sig")
    qe.set_index('set_id', inplace = True)

    for set_id, row in qe.iterrows():
        id_directory = os.path.join(homeDirectory, set_id)
        createDirectory(log,id_directory, f"# INFO - Directory for paramater set [{set_id}] created at [{id_directory}]", False)
    
        #Make sure there is a directory to process
        input_files = os.path.join(id_directory, "input_files")
        createDirectory(log,input_files, "# WARN - No directory found for input files.", True)
        directories = sorted([directory for directory in os.listdir(input_files) if os.path.isdir(os.path.join(input_files, directory)) and not directory.startswith(".") and not os.path.isfile(os.path.join(os.path.join(input_files, directory), "INCOMPLETE.txt"))])
    
        for refcode in directories:
            setid_refcode = set_id + " " + refcode
            if setid_refcode in df['[REFCODE]'].values:
                printToLog(log,"# INFO - Compound ["+ setid_refcode +"] Already present in sheet")
            else:           
                printToLog(log,"# INFO - Compound ["+ setid_refcode +"] Appending to sheet")
                df = pd.concat([df, pd.DataFrame({"[REFCODE]": [setid_refcode]})], ignore_index=True)
                processedCount += 1

    df.to_csv(localSheet, index=False)
    printToLog(log,"# INFO - Successfully appended ["+str(processedCount)+"] compounds to sheet at ["+str(localSheet)+"]")

# Extracts data from summary files and updates a local .csv
def updateCSV(log):
    printToLog(log,"# INFO - Attempting to update sheet at ["+str(localSheet)+"]")
    if not os.path.isfile(localSheet):
        printToLog(log, "# WARN - No local .csv found.")
        quit()

    if not os.path.isfile(qe_params):
        printToLog(log,"# WARN - No qe_params.csv found.")
        quit()

    qe = pd.read_csv(qe_params, encoding="utf-8-sig")
    qe.set_index('set_id', inplace = True)
    for set_id, row in qe.iterrows():
        processedCount = 0
        printToLog(log,f"# INFO - Processing set_id [{set_id}]")
    
        id_directory = os.path.join(homeDirectory, set_id)
        createDirectory(log,id_directory, f"# INFO - Directory for paramater set [{set_id}] created at [{id_directory}]", False)
    
        summary_files = os.path.join(id_directory, "summary_files")
        createDirectory(log,summary_files, "# INFO - No directory found for summary files, creating at ["+str(summary_files)+"]", False)
        summaries = sorted([file for file in os.listdir(summary_files) if file.endswith('_summary.txt') and os.path.isfile(os.path.join(summary_files, file))])
        printToLog(log,"# INFO - The following summaries are available ["+str(summaries)+"]")
    
        df = pd.read_csv(localSheet)  
        df.set_index('[REFCODE]', inplace = True)
        df = df.astype(str)
        
        for summary in summaries:
            refcode = os.path.splitext(summary)[0].replace("_summary", "")
            setid_refcode = set_id + " " + refcode
            printToLog(log,"# INFO - Compound ["+ setid_refcode +"] Processing output data")
    
            with open(os.path.join(summary_files, summary), "r") as file:            
                read = file.read()
                lines = read.splitlines()
                
                for line in lines:
                    if not len(line) == 0 and not line.startswith("#") and not line.startswith("_"):
                        value = line[line.find("=")+1:].strip()
                        name = line[:line.find("=")-1].strip()
                        writeCSV(df, setid_refcode, "["+str(name)+"]", str(value))
                processedCount += 1
        
        df = df.replace("nan", "")
        df.to_csv(localSheet) # Update local csv
        printToLog(log,"# INFO - Successfully updated data in sheet at ["+str(localSheet)+"] for ["+str(processedCount)+"] compounds")

# References and updates a local .csv to submit requests to slurm, only runs calculations not already flagged as batched
# Batches up to 'batchCount' (16) every run to avoid requesting too many resources at once 
def batchCalculations(log, batchCount):
    printToLog(log,"# INFO - Attempting to batch ["+str(batchCount)+"] calculations")
    if not os.path.isfile(localSheet):
        printToLog(log, "# WARN - No local .csv found.")
        quit()

    if not os.path.isfile(qe_params):
        printToLog(log,"# WARN - No qe_params.csv found.")
        quit()

    processedCount = 0
    qe = pd.read_csv(qe_params, encoding="utf-8-sig")
    qe.set_index('set_id', inplace = True)
    for set_id, row in qe.iterrows():
        if processedCount >= batchCount:
            break
        
        printToLog(log,f"# INFO - Processing set_id [{set_id}]")
    
        id_directory = os.path.join(homeDirectory, set_id)
        createDirectory(log,id_directory, f"# INFO - Directory for paramater set [{set_id}] created at [{id_directory}]", False)
    
        #Make sure there is a directory to process
        input_files = os.path.join(id_directory, "input_files")
        createDirectory(log,input_files, "# WARN - No directory found for input files.", True)
        directories = sorted([directory for directory in os.listdir(input_files) if os.path.isdir(os.path.join(input_files, directory)) and not directory.startswith(".") and not os.path.isfile(os.path.join(os.path.join(input_files, directory), "INCOMPLETE.txt"))])
    
        if len(directories) == 0:
            printToLog(log,"# WARN - No directories found in ["+ input_files + "]")
            quit()
        else:
            printToLog(log,f"# INFO - [{len(directories)}] directories found at [{input_files}]")

        df = pd.read_csv(localSheet)  
        df.set_index('[REFCODE]', inplace = True)
  
        for refcode in directories:
            setid_refcode = set_id + " " + refcode
            if processedCount < batchCount:
                printToLog(log,"# INFO - Processing compound with refcode ["+ setid_refcode +"]")
                if not str(df.at[setid_refcode, "[BATCH_started]"]) == "True":
                    printToLog(log,"# INFO - Compound with refcode ["+ setid_refcode +"] not previously run, attempting to batch")
                    refcodeDirectory = os.path.join(input_files, refcode)
                    
                    QE_SUB = os.path.join(refcodeDirectory, "QE_SUB")
                    batch_path = os.path.join(refcodeDirectory, refcode+"_batch.txt")
    
                    batchCommand = f"module load {param_modules}; cd {refcodeDirectory}; sbatch QE_SUB"
                    if os.path.exists(QE_SUB):
                        try:
                            now = str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                            subprocess.call(batchCommand,shell=True)
                            
                            with open(batch_path, "a") as batch:
                                writeCSV(df, setid_refcode, "[BATCH_started]", True)
                                writeCSV(df, setid_refcode, "[BATCH_start_time]", now)
                                writeCSV(df, setid_refcode, "[BATCH_location]", param_location)
                                
                                print("\n# -Batch data\n", file=batch)
                                print("BATCH_started = "+str(True), file=batch)
                                print("BATCH_start_time = "+str(now), file=batch)
                                print("BATCH_location = "+str(param_location), file=batch)
                                
                                printToLog(log,"# INFO - Successfully batched calculation for compound ["+setid_refcode+"] at ["+str(now)+"] on ["+str(param_location)+"]")
                                processedCount += 1
                        except subprocess.CalledProcessError as e:
                            printToLog(log,"# WARN - Error batching calculation for compound with refcode ["+setid_refcode+"]")
                            printToLog(log,str(e))
                    else:
                        printToLog(log,"# WARN - QE_SUB not present for compound with refcode ["+setid_refcode+"]")
                else:
                    printToLog(log,"# INFO - Compound with refcode ["+setid_refcode+"] has been previously batched at ["+str(df.at[setid_refcode, "[BATCH_start_time]"])+"] on ["+str(df.at[setid_refcode, "[BATCH_location]"])+"]")  
            else:
                break
        df.to_csv(localSheet)#Update local csv
        printToLog(log,"# INFO - ["+str(processedCount)+"] Calculations successfully batched.")
    
    if processedCount < batchCount:
        printToLog(log,"# INFO - No more calculations to batch!")
    getQueueLength(log)
    

        

















# Unused


#def initSheet(log):
#    if os.path.isfile(localSheet):
#        printToLog(log,"# INFO - Removing existing local file ["+ localSheet + "]")
#        os.remove(localSheet)# Clear current local copy
        
#    printToLog(log,"# INFO - Creating new sheet ["+ localSheet + "]")
#    with open(localSheet, 'a') as file:
#        file.write("[REFCODE]")
#    df = pd.read_csv(localSheet)
#    df = pd.concat([df, pd.DataFrame({"[REFCODE]": ["init"]})], ignore_index=True)    

#    df["[BATCH_location]"] = ["Abyss"]    
#    df["[BATCH_started"] = ["True"]
#    df.to_csv(localSheet, index=False)

# Unused - too slow to be helpful
#def isolateFailuresInCSV(log):
#    failureCount = 0    
#    printToLog(log,"# INFO - Attempting to find failures in sheet at ["+str(localSheet)+"]")

#    if not os.path.isfile(localSheet):
#        printToLog(log, "# WARN - No local .csv found.")
#        quit()

#    failures = os.path.join(homeDirectory, "failures")
#    createDirectory(log, failures, "# WARN - No directory found for failed input files.", False)
#    previous = [directory for directory in os.listdir(failures) if os.path.isdir(os.path.join(failures, directory)) and not directory.startswith(".")]
#    
#    input_files = os.path.join(homeDirectory, "input_files")
#    createDirectory(log, input_files, "# WARN - No directory found for input files.", True)
#    directories = [directory for directory in os.listdir(input_files) if os.path.isdir(os.path.join(input_files, directory)) and not directory.startswith(".")]

#    directories = sorted(directories)
#    printToLog(log,"# INFO - The following input directories are available ["+str(directories)+"]")
    
#    df = pd.read_csv(localSheet)
#    df.set_index('[REFCODE]', inplace = True)
#    for refcode, row in df.iterrows():        
#        if previous.__contains__(refcode):
#            continue
#        if not str(row["[BATCH_started]"]) == "True":
#            continue
        
#        failed = False
#        refcodeDirectory = os.path.join(input_files, refcode)
        
#        if directories.__contains__(refcode):            
#            if not str(row["[BATCH_done]"]) == "True":
#                if not isQueued(log, refcode):
#                    printToLog(log,"# INFO - Compound ["+ refcode +"] Did not finish")
#                    failed = True
#            else:
#                if not str(row["[PWSCF_done]"]) == "True":
#                    printToLog(log,"# INFO - Compound ["+ refcode +"] PWSCF Did not finish")
#                    failed = True

#                if not str(row["[GIPAW_done]"]) == "True":
#                    printToLog(log,"# INFO - Compound ["+ refcode +"] GIPAW Did not finish")
#                    failed = True

#                if str(row["[PWSCF_finalEnergy]"]) == "nan":
#                    printToLog(log,"# INFO - Compound ["+ refcode +"] Did not converge")
#                    failed = True
#        else:
#            printToLog(log,"# INFO - Compound ["+ refcode +"] Not present on this cluster")
#        if failed:
#            shutil.copytree(refcodeDirectory, os.path.join(failures, refcode))
#            failureCount += 1
        
#    printToLog(log,"# INFO - Determined that ["+str(failureCount)+"] calculations have failed")

# Try to determine location from the host name - not always very successful - unused
# "local_log" Should in theory be the log of whatever script is currently trying to use getLocation() 

#local_log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
#utils = os.path.join(homeDirectory, "utils")
#location = os.path.join(utils, "location.txt")
#if not os.path.exists(location): 
#    with open(location, "a") as file:
#        printToLog(log,local_log," --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - No location data found. Attempting to retreive.")    
#        try:
#            out = subprocess.check_output(['hostname'],shell=True)
#            out = out.decode("utf-8").strip()
#    
#            print(out,file=file)
#            printToLog(log,local_log,"# INFO - Location determined to be ["+str(out)+"], saved to ["+str(location)+"]")
#            printToLog(log,local_log,"# INFO - Override manually by changing the contents of ['location.txt']")
#        except subprocess.CalledProcessError as e:
#            printToLog(log,local_log,"# INFO - Error retreiving location data.")
#            printToLog(log,local_log,str(e))

#def getLocation():
#    with open(location, "r") as file:
#        return file.read().strip() # Simply reads the file


