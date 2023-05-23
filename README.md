## Podacris Genome Annotation (PGA)

TODO: Link manuscript when published
### Genomic signatures of selection underpinning extremely rapid adaptation in island lizards
#### Authors: Maria Novosolov, Anthony Herrel, Anamaria Stambuk, Jesper Stenderup, James S. Santangelo, Iva Sabolić, Kristian E. Hanghøj, Thorfinn S. Korneliussen, Debora Y. C. Brandt, Rasmus Heller, Rasmus Nielsen, Morten E. Allentoft



### Description

This repository contains the code used for annotation scaffold-level assembly of the Italian Wall Lizard
(*Podacris siculus*). The pipeline starts with input of the genome FASTA file. It outputs the final
functional annotation in GFF3 format

In addition to the genome FASTA, the following dependencies/resources need to be obtained prior to running the pipeline.

1. The pipeline uses GeneMark, which is licensed softwared. The license key can be obtained from 
    [here](http://topaz.gatech.edu/Genemark/license_download.cgi).  
    1. Include the `.gm_key` in the [workflow/](./workflow) directory containing the main Snakefile used for running the pipeline
2. A copy of the RepBase repeat library (I used `RepBaseRepeatMaskerEdition-20181026.tar.gz`) from [here](https://www.girinst.org/repbase/). 
    A license is required to use this database so it could not be included. Once obtained, it needs to be included in [resources/](./resources) directory
3. A copy of the `InterProScan` data needs to be obtained and included in the [resources/](./resources) directory. 
    The following commands can be run to download and setup the database:

    ```
    wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.61-93.0/alt/interproscan-data-5.61-93.0.tar.gz
    tar -xvzfp interproscan-data-5.61-93.0.tar.gz
    chmod -R 777 interproscan-data-5.61-93.0
    ```

#### Using the pipeline

This pipeline requires `Conda` and `Singularity`:

A minimal installation of `Conda` (i.e., `Miniconda`) can be installed by following the instructions for your platform [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

Installation of `Singularity` requires Admin privileges, but using `Singularity` to run pre-created containers does not. 
Installation instructions can be found [here](https://docs.sylabs.io/guides/latest/admin-guide/). 
All `Singularity` containers used in this pipeline are avalaible in [this public repository](https://cloud.sylabs.io/library/james-s-santangelo), 
though they will be automatically pulled and executed by the pipeline.
Assuming `Conda` is installed, the this repository's Conda environment can be replicated by running the following command:

```
conda env create -f environment.yaml -n prg
```

This will create a `Conda` environment named prg containing a minimal set of dependencies required to run the pipeline. 

After activating the environment (`conda activate dcg`), the pipeline can be executed from the workflow directory by running a command that looks something like:

```
snakemake --use-conda --use-singularity --singularity-args "--bind <path> --bind <path/to/interproscan/data>:/opt/interproscan-5.61-93.0/data" --configfile ../config/<configfile> --notemp -j <cores>
```

for local execution. Here, `<path>` is the path on the cluster from which files will be read/written (e.g., `/scratch`), `<path/to/interproscan/data>` is the full path to the interproscan data in [resources/](./resources), `<configfile>` is one of the configfiles in the config directory that needs to be modified to match the paths on your system, and `<cores>` is the number of cores available for executing parallel processes.

