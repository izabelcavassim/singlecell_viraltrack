# Singlecell viraltrack analyses

The present directory exemplify how to detect Viral DNA in single cell human data using the software [Viral-track](https://github.com/PierreBSC/Viral-Track).

To run the analyses above you will require to have conda installed and will also be required to use the conda environment (single_cell.yaml) shared in this repository [here](https://github.com/izabelcavassim/singlecell_viraltrack/blob/main/single_cell.yaml).

To install the conda environment from the yaml file you can type:
```
conda env create -f single_cell.yaml

# then activate your environment
conda activate single_cell
```
Once we have it installed, you can now dig into the workflow (scripts_singlec_UMI.py), found [here](https://github.com/izabelcavassim/singlecell_viraltrack/blob/main/scripts_singlec_UMI.py) written with the software [gwf](https://gwf.app/).

It is expected that the files are in specific directories. You would want to modify the workflow paths to assure that you have the same organization of files as the workflow requires. Once this is all organized then you can run the pipeline by typing:

```
 gwf -f scripts_singlec_UMI.py run
```
These will submit the jobs to the cluster. The cluster uses SGE engine, so bare in mind that. 

To see the status of your job you can type:
```
 gwf -f scripts_singlec_UMI.py status
```

I have written a tutorial on how to use and install gwf in the hoffman2 cluster (at UCLA). You can find it [here](https://github.com/izabelcavassim/gwf_UCLA_cluster/blob/master/gwf_tutorial.md).

Feel free to reach out if any problems arise (izabelcavassim    at.  gmail)


