#Generic utility functions used by all scripts

import csv
import datetime
import time
import os
import shutil
import subprocess
import numpy as np

logs = os.path.join(os.getcwd(), "logs")
if not os.path.exists(logs):
    os.makedirs(logs)

#Prints and logs in one, convention I personally like
def printToLog(log, info):
    info = str(info)
    time = ""
    if not info.startswith(" ---"):
        time = str(datetime.datetime.now().strftime("[%H:%M:%S] "))
    print(time+str(info))
    with open(os.path.join(logs, log), "a") as log:
        log.write(time+str(info) + "\n")
        
#Create directory if it doesn't exist. Optionally crash deliberately if doesn't exist
def createDirectory(log, path, text, exit):
    if not os.path.exists(path):
        printToLog(log, text + " ["+ path + "]")
        os.makedirs(path)
        if exit:
            quit()

#Remove directory is it exists
def removeDirectory(log, path, text):
    if os.path.exists(path):
        printToLog(log, text + " ["+ path + "]")
        shutil.rmtree(path)

#Write an entry to the local csv
def writeCSV(df, refcode, location, value):
    if not value == "":
        df.loc[refcode, location] = value

def cellVolume(cell_a, cell_b, cell_c, cell_α, cell_β, cell_γ):
    return cell_a * cell_b * cell_c * np.sqrt((1 - (np.cos(cell_α) ** 2) - (np.cos(cell_β) ** 2) - (np.cos(cell_γ)) ** 2) + (2 * np.cos(cell_α) * np.cos(cell_β) * np.cos(cell_γ)))

#Get and decode the current slurm queue. Can be read like a file
def getQueue(log):
    printToLog(log,"# INFO - Attempting to retrieve current slurm queue.")
    try:
        out = subprocess.check_output(['squeue --me'],shell=True)
        out = out.decode("utf-8")
        return out
    except subprocess.CalledProcessError as e:
        printToLog(log,"# INFO - Error retreiving slurm queue.")
        printToLog(log,str(e))

#Get the length of the current slurm queue
def getQueueLength(log):
    printToLog(log,"# INFO - Attempting to get the length of current slurm queue.")

    length = 0
    lines = getQueue(log).splitlines()
    for line in lines:
        if "_SUB" in line:
            length += 1
        printToLog(log, line)
    printToLog(log,"# INFO - Slurm queue contains ["+str(length)+"] batched calculations.")
    return length

#Check if a specific refcode is in the queue
def isQueued(log, refcode):
    printToLog(log,"# INFO - Compound ["+refcode+"] Checking queue")
    lines = getQueue(log).splitlines()
    for line in lines:
        if refcode+"_SUB" in line:
            printToLog(log,"# INFO - Compound ["+refcode+"] is currently queued.")
            return True
    printToLog(log,"# INFO - Compound ["+refcode+"] is not currently queued.")
    return False