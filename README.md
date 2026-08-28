Designed for use on Digital Research Alliance of Canada clusters.

# To activate:
Run the following commands:
```
git clone https://github.com/Jaydee8652/EspressoBatch.git
```
```
cd EspressoBatch/
```
```
module load scipy-stack/2023b
```
```
pip install -r requirements.txt
```


Edit 'utils/params.py' 


### User Specific Parameters

-Replace the value for "param_email" with the name of an email you control. Slurm events will be forwarded to this email.

-Replace the value for "param_account" with your account on the cluster

-The value for "param_slurmVerbosity" can be any supported slurm --mail-type (https://slurm.schedmd.com/sbatch.html)

### Cluster Specific Parameters (See node characteristics for your cluster: https://docs.alliancecan.ca/wiki/Rorqual/en)

-Replace the value for "param_location" with the human readable name of the cluster (ie "Rorqual")

-Replace the value for "param_cores" with the number of cores on each node 

-Replace the value for "param_memory" with the available memory of each node in G



### You can now run calculations!
It should be noted that the scipy-stack must to be reloaded every time you connect to the cluster. Creating an alias to do this automatically is recommended.

```
module load scipy-stack/2023b
```

# Git integration
By default, batch status and calculation outputs are saved to a .csv locally. Optionally, github integration can be enabled. This data will instead be saved to a defined repository, allowing the same global database to be referenced across multiple clusters. 

To activate this feature, enter 'utils/params.py' 


-Replace the value for "param_repo" with the name of a github repository you control
-Replace the value for "param_token" with an auth token that has read and write permissions on said repo

In the repository:

-Create a new .csv file called 'sheet.csv' in the home directory
-Create a new .txt file called 'sheet_flag.txt' in the home directory, containing the string 'True'



# Usage instructions:



## cif_sort.py
```
python3 cif_sort.py
```
Intended to quickly filter a directory of .cif files by several characteristics listed below. 

If the provided .cif files are from the Cambridge Structural Database, a .csv of additional characteristics can be obtained by saving a selection of structures as a TAB separated values table in Conquest. This can then be converted with excel or another .csv manager. This file can be provided at 'EspressoBatch/structure_data.csv' to allow for more robust filtering. Without a provided 'structure_data.csv', filtering settings that require it will be forcibly disabled if selected. 

On the first run this script will create the 'original_cifs' directory. Place .cif files in 'original_cifs' and rerun, the user will then be presented with different filtering settings. Any combination of these filters can be run through an integer input. 

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

Once the sort is complete .cif files are moved to directories within the 'cifs' directory corresponding to the reason for discard, 'original_cifs' is also moved to this directory. Running the script again will create a new 'original_cifs' to start the process again. Previously sorted .cif files will not be overriden, allowing for a new directory of .cif files to be added to the existing filtered dataset.



## sanity_check.py
```
python3 sanity_check.py
```
This script is intended for use in crystal structure prediction, allowing the energy of generated .cif files to be quickly assessed and unreasonable structures discarded. Streamlines the process of performing "sanity check" calculations (scf calculations with a single K point to estimate energy).

Not considered necessary for .cifs obtained from the Cambridge Structural Database, all CSD structures will likely have reasonable energies. 

Creates the 'Sanity_Input_Files' directory and presents the user with 4 processes to run. Any of these processes can be run with an integer input, but they should ideally be performed in order.

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


## qe_cif2cell.py
```
python3 qe_cif2cell.py
```
Will generate quantum-espresso input files with calculation parameters defined in 'qe_params.csv' from .cifs in 'cifs/validated'. Automatically runs test calculations and creates slurm batch files based on the the projected resource use.

On the first run this script will create 'qe_params.csv' with default settings. Adding new rows to 'qe_params.csv' will create multiple sets of calculations, each set will be given its own directory. Input files will be created for all sets with that sets defined parameters.

| set_id | test_time | ecutwfc | ecutrho_factor | conv_thr | q_gipaw | calculation | volume_cap | atoms_to_optimise |
| ------ | --------- | ------- | -------------- | -------- | ------- | ----------- | ---------- | ----------------- |
| MAIN   | 10        | 55.0    | 8.0            | 1.D-6    | 0.01    | relax       | 0          | H                 |


 - set_id (default: MAIN)
 Human readable name for this "set". 
 Used as the name of the directory and will be appended to .cif refcodes in the slurm queue and output spreadsheet.

 - test_time (default: 10)
Time taken for a test calculation to complete (in seconds). 
Used to estimate the amount of time creating a series of input files will take. Slightly varies with calculation difficulty.

 - ecutwfc (default: 55.0)
*Default quantumespresso parameter

 - ecutrho_factor: (default: 8.0)
Multiplier of ecutwfc to produce quantumespresso parameter ecutrho.

 - conv_thr (default: 1.D-6)
*Default quantumespresso parameter

 - q_gipaw (default: 0.01)
*Default quantumespresso parameter

 - calculation (default: relax)
Type of calculation to perform, see quantumespresso manual. 
Only relax is fully supported by post processing scripts.

 - volume_cap (default: 0)
The minimum volume at which calculations should be given more time to process.
Above this number, the time requested from slurm increased from 1 day to 3 days.

 - atoms_to_optimise: (default: H)
List of atom types to optimise, all other atom types will be frozen.
Should be seperated by " " (ie: H C O N Co).
Enter "None" to optimise no atom positions and "All" to optimise all atom positions.
 

Presents the user with 2 processes to run. Only one of these processes can be run at once, through an integer input.

 - 1:
Run test calculations in a slurm job array

Determines appropriate groupings of input files based on the time each test calculation will take and submits them to slurm. Empty directories will be created first, and the job array will populate them.

 - 2: 
Run test calculations in current session

Running test calculations in the current session is more brittle, as if the current session ends the process will stop. It is useful if only a few dozen input files are required as there is no wait for slurm resource allocation.
 




## batch_control.py
```
python3 batch_control.py
```
Displays the current slurm queue, determines the number of previously batched calculations, and then presents the user with 3 processes to run. Any combination of these processes can be run through an integer input.

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
 
 
 It should be noted that if 2 and 3 attempt to modify a refcode not in the .csv they may crash.
