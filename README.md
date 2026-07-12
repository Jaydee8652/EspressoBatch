Designed for use on Digital Research Alliance of Canada clusters.

To activate, run the following commands:
```
$ git clone https://github.com/Jaydee8652/EspressoBatch.git
$ cd EspressoBatch/
$ pip install -r requirements.txt
$ module load scipy-stack
```

By default, batch status and calculation outputs are saved to a .csv locally. Optionally, github integration can be enabled. This data will instead be saved to a defined repository, allowing the same global database to be referenced across multiple clusters. Activating github integration:

In 'utils.git_utils'
```
-Replace the value for "REPO" with the name of a github repository you control
-Replace the value for "TOKEN" with an auth token that has read and write permissions on said repo
```
In the repository:
```
-Create a new .csv file called 'sheet.csv' in the home directory
-Create a new .txt file called 'sheet_flag.txt' in the home directory, containing the string 'True'
```

Finally run the following command to populate the .csv with headers:
```
$ python3 _init.py
```



Usage instructions:

```
$ python3 cif_sort.py
```
On first usage will create the 'Original_CIFs' directory and request a structure_data.csv. The structure_data.csv can be obtained from the CSD, by saving a selection of structures as a TAB seperated values table and converting with excel or other .csv manager. Place corresponding .cif files in 'Original_CIFs', also obtained from the CSD and rerun. The sort should complete within a minute.

It should be noted that .cifs with "structure-formula mismatch" may contain an unreported cocrystal/solvent. Manual inspection of these files is encouraged.

```
$ python3 qe_cif2cell.py
```
Will generate quantum-espresso input files, automatically run a test calculation, and then update the slurm request according to the projected resource use.

```
$ python3 batch_control.py
```
On first usage will attempt to determine the name of the local cluster, saving it to 'location.txt' to be referenced by other scripts. This string can be changed manually if necessary. 

Retreives the current slurm queue, determines the number of previously batched calculations, and then presents the user with 3 processes to run. Any combination of these processes can be run through an integer input.
```
 1
 Append the refcode of all local input directories to a .csv stored locally / on github

 2 
 Extract data from local summary files and updates a .csv stored locally / on github
 Intended to be run after a series of calculations have finished, inclusion in the workflow here allows the previous 
 batch to be processed when a new one is requested. Summary files are produced at the end of an sbatch calculation by 'post_processing.py'

 3 
 References and updates a .csv stored locally / on github to submit requests to slurm. Will only run calculations not flagged as previously batched.
 Batches 'batchCount' every run to avoid requesting too many resources at once. 'batchCount' by default is the number of calculations that would lead to a slurm queue length of 16. 

 0
 "Speed dial" for all processes in sequence
 ```
 
 It should be noted that if 2 and 3 attempt to modify a refcode not in the .csv they may crash, user discression is advised.
