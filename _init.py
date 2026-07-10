# Run to replace the current sheet stored on github (if any) with a new blank sheet

#Imports
import os
import csv
import re
import shutil
import sys
import pandas as pd
import datetime
import time
from utils.generic_utils import printToLog as pl, createDirectory as cd, removeDirectory as rd
from utils.git_utils import downloadCSV, uploadCSV, getLocation, verify, initSheet

#Functions
def printToLog(info):#Prints and logs in one, convention I personally like
    pl(log, info)
def createDirectory(path, text, exit):
    cd(log, path, text, exit)
def removeDirectory(path, text):
    rd(log, path, text)

#Main
log = str(os.path.basename(sys.argv[0]).split(".")[0]+".log")
homeDirectory = os.getcwd()#Directory where we are
printToLog(" --- \n"+str(datetime.datetime.now().strftime("[%H:%M:%S] "))+"# INFO - Starting new "+str(os.path.basename(sys.argv[0]).split(".")[0])+" process in ["+ homeDirectory + "]")    

initSheet(log)
uploadCSV(log)

