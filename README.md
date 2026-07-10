Designed for use on Digital Research Alliance of Canada clusters.

To activate, run the following commands:

```
$ git clone https://github.com/Jaydee8652/EspressoBatch.git
$ cd EspressoBatch/
$ pip install -r requirements.txt
$ module load scipy-stack
```

By default, batch status and calculation outputs are saved to a .csv locally. Optionally, github inegration can be enabled to instead save this data to a repository, allowing the same gloabl database to be referenced from multiple clusters.

Activating github integration:

In utils.git_utils:
```
-Replace "REPO" with the name of a github repository you control
-Replace "TOKEN" with an auth token that has read and write permissions on said repo
```

Use:

Run cif_sort.py first, this will create the Original_CIFs file and request structure_data.csv, both can be obained from the CSD.

Run qe_cif2cell.py to generate quantum-espresso input files.
