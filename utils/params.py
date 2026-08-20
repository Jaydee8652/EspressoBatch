# Params - can be changed

#MODULE VERSIONS
param_modules = "StdEnv/2023 quantumespresso/7.3.1 scipy-stack/2023b xtb/6.6.1 openbabel/3.1.1"

#SLURM
param_email = "EMAIL@gmail.com"
param_account = "def-ACCOUNT"
param_slurmVerbosity = "ALL"

#CLUSTER SPECIFIC
param_location = "Rorqual" # String name for cluster
param_cores = 192 # See node characteristics
param_memory = 750 #in G - See node characteristics

#GITHUB
param_repo = "REPO_NAME" #https://github.com/Jaydee8652/REPO_NAME
param_token = "github_pat_0000000000000000000000000000000000000000000000000000000000000000000000000000000000" #Must have permissions on repo
param_sheetPath = 'sheet.csv'
param_flagPath = 'sheet_flag.txt'