#Generic utility functions used by all scripts

import csv
import datetime
import time
import os
import shutil

# Params - can be changed
modules = "StdEnv/2023 quantumespresso/7.3.1 scipy-stack/2023b xtb/6.6.1"


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

def getModules():
    return modules 

