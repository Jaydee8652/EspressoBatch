#Generic utility functions used by all scripts
import csv
import datetime
import time
import os
import shutil
import subprocess
import numpy as np
import re

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