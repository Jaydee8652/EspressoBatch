Designed for use on Digital Research Alliance of Canada clusters.

# To activate:

Run the following commands:
```
$ git clone https://github.com/Jaydee8652/EspressoBatch.git
$ cd EspressoBatch/
$ pip install -r requirements.txt
$ module load scipy-stack/2023b
```

Once downloaded, enter 'utils.params.py' 

```
-Replace the value for "param_email" with the name of an email you control. Slurm events will be forwarded to this email.
-Replace the value for "param_account" with your account on the cluster
-The value for "param_slurmVerbosity" can be any supported slurm --mail-type (https://slurm.schedmd.com/sbatch.html)
```

# You can now run calculations!

By default, batch status and calculation outputs are saved to a .csv locally. Optionally, github integration can be enabled. This data will instead be saved to a defined repository, allowing the same global database to be referenced across multiple clusters. 

To activate this feature, enter 'utils.params.py' 

```
-Replace the value for "param_repo" with the name of a github repository you control
-Replace the value for "param_token" with an auth token that has read and write permissions on said repo
```
In the repository:
```
-Create a new .csv file called 'sheet.csv' in the home directory
-Create a new .txt file called 'sheet_flag.txt' in the home directory, containing the string 'True'
```


# Usage instructions:

```
$ python3 cif_sort.py
```
On first run will create the 'Original_CIFs' directory, place .cif files in 'Original_CIFs' and rerun. The user will be presented with different filtering settings, some will be disabled without a 'structure_data.csv'. Any combination of these filters can be run through an integer input. 

If the .cifs are from the Cambridge Structural Database, the 'structure_data.csv' can be obtained from the CSD. This can be done by saving a selection of structures as a TAB separated values table and converting with excel or other .csv manager. Without a provided 'structure_data.csv' the filtering settings that require it will be forcibly disabled if selected.

```
$ python3 qe_cif2cell.py
```
Will generate quantum-espresso input files, automatically run a test calculation, and then update the slurm request according to the projected resource use.

```
$ python3 batch_control.py
```
On first run will attempt to determine the name of the local cluster, saving it to 'location.txt' to be referenced by other scripts. This string can be changed manually if necessary. 

Displays the current slurm queue, determines the number of previously batched calculations, and then presents the user with 3 processes to run. Any combination of these processes can be run through an integer input.
```
 - 1:
Append the refcode of all local input directories to a .csv stored locally / on github

Input directories are produced by qe_cif2cell.py

 - 2: 
Extract data from local summary files and update a .csv stored locally / on github
Intended to be run after a series of calculations have finished, inclusion in the workflow here allows the previous 
batch to be processed when a new one is requested.

Summary files are produced at the end of an sbatch calculation by 'post_processing.py'

 - 3: 
References and updates a .csv stored locally / on github to submit requests to slurm.
Creates 'REFCODE_batch.txt' to store the time and location of the batch.
Will only run calculations not flagged as previously batched.

Batches 'batchCount' every run to avoid requesting too many resources at once.
'batchCount' by default is the number of calculations that would lead to a slurm queue length of 16. 

Displays the final slurm queue once batching is complete

 - 0:
"Speed dial" for all processes in sequence
 ```
 
 It should be noted that if 2 and 3 attempt to modify a refcode not in the .csv they may crash.
