Designed for use on Digital Research Alliance of Canada clusters.

# To activate:
Run the following commands:
```
$ git clone https://github.com/Jaydee8652/EspressoBatch.git
$ cd EspressoBatch/
$ pip install -r requirements.txt
$ module load scipy-stack/2023b
$ cd utils
```

Edit 'params.py' 

```
-Replace the value for "param_email" with the name of an email you control. Slurm events will be forwarded to this email.
-Replace the value for "param_account" with your account on the cluster
-The value for "param_slurmVerbosity" can be any supported slurm --mail-type (https://slurm.schedmd.com/sbatch.html)
```

```
$ cd ..
```

On the first run of some scripts they will attempt to determine the name of the local cluster, saving it to 'utils/location.txt'. This string can be changed manually if necessary. 

### You can now run calculations!
By default, batch status and calculation outputs are saved to a .csv locally. Optionally, github integration can be enabled. This data will instead be saved to a defined repository, allowing the same global database to be referenced across multiple clusters. 

To activate this feature, enter 'utils/params.py' 

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

## cif_sort.py
```
$ python3 cif_sort.py
```
Intended to quickly filter a directory of .cif files by several characteristics listed below. 

If the provided .cif files are from the Cambridge Structural Database, a .csv of additional characteristics can be obtained by saving a selection of structures as a TAB separated values table in Conquest. This can then be converted with excel or another .csv manager. This file can be provided at 'EspressoBatch/structure_data.csv' to allow for more robust filtering. Without a provided 'structure_data.csv', filtering settings that require it will be forcibly disabled if selected. 

On the first run this script will create the 'original_cifs' directory. Place .cif files in 'original_cifs' and rerun, the user will then be presented with different filtering settings. Any combination of these filters can be run through an integer input. 
```
 - 1:
Discard structures with r factor greater than [{rCap} (default: 10)]
*(This data is stored in 'structure_data.csv')

 - 2: 
Discard structures flagged as disordered by their CSD author
*(This data is stored in 'structure_data.csv')

 - 3: 
Discard structures with volume greater than [{volumeCap} (default: 6000)]

 - 4:
Discard structures without structural data

 - 5:
Discard structures without hydrogen data

 - 6:
Discard structures with incomplete hydrogen data

 - 7:
Discard structures with unreported cocrystals/solvent

 - 0:
"Speed dial" for all filters in sequence
 ```

Once the sort is complete .cif files are moved to directories within the 'cifs' directory corresponding to the reason for discard, 'original_cifs' is also moved to this directory. Running the script again will create a new 'original_cifs' to start the process again. Previously sorted .cif files will not be overriden, allowing for a new directory of .cif files to be added to the existing filtered dataset.



## sanity_check.py
```
$ python3 sanity_check.py
```
This script is intended for use in crystal structure prediction, allowing the energy of generated .cif files to be quickly assessed and unreasonable structures discarded. Streamlines the process of performing "sanity check" calculations (scf calculations with a single K point to estimate energy).

Not considered necessary for .cifs obtained from the Cambridge Structural Database, all CSD structures will likely have reasonable energies. 

Creates the 'Sanity_Input_Files' directory and presents the user with 4 processes to run. Any of these processes can be run with an integer input, but they should ideally be performed in order.
```
 - 1:
Generates quantum-espresso sfc input files from .cifs in 'cifs/validated' to perform a "sanity check".

 - 2: (MUST BE RUN ON HEAD NODE)
Batches all sanity check calculations to slurm. Calculations are grouped into jobs dynamically, with an initial target of 1000 calculations per job submitted.
These calculations will run sequentually within a job. If calculations are particularly demanding, or speed is a concern, the grouping can reduced to split
the task across more slurm jobs.

All calculations will be batched at once, regardless of group size as long as the number of jobs requested does not exceed 16.
Automatcially logs which calculations have been batched in 'Sanity_Input_Files/sanity_sheet.csv' and will not repeat work.

Summary files are produced at the end of each calculation by 'extract_energy.py'
The final energy is saved to 'Sanity_Input_Files/sanity_sheet.csv' in Ry and kJ mol⁻¹ molecule⁻¹. 

 - 3: 
Should only be run after all calculations are complete. Determines the relative energy of all outputs in kJ mol⁻¹ molecule⁻¹.
The lowest energy output is used to 'zero' all other energies. 

The energy of each refcode can then be inspected manually, structures with reasonable energies can be marked to be retained.
This is done by entering 'True' in the [validated] column for a given row in 'Sanity_Input_Files/sanity_sheet.csv'.

 - 4: 
Discards .cifs from 'cifs/validated' not marked in 'Sanity_Input_Files/sanity_sheet.csv' with [validated] = 'True'. Discarded structures are
saved to 'cifs/high_energy' and a backup of 'cifs/validated' is created so that the test can be run again with different parameters if desired.
 ```



## qe_cif2cell.py
```
$ python3 qe_cif2cell.py
```
Will generate quantum-espresso input files from .cifs in 'cifs/validated', automatically run a test calculation, and create a batch file according to the projected resource use.
Takes a list user input of atom types to optimise, all other atom types will be frozen.



## batch_control.py
```
$ python3 batch_control.py
```
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

Batches [batchCount] every run to avoid requesting too many resources at once.
[batchCount] by default is the number of calculations that would lead to a slurm queue length of 16. 

Displays the final slurm queue once batching is complete

 - 0:
"Speed dial" for all processes in sequence
 ```
 
 It should be noted that if 2 and 3 attempt to modify a refcode not in the .csv they may crash.
