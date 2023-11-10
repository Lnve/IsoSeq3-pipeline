# IsoSeq3-pipeline
The pipeline provides an automated Snakemake workflow for the identification and annotation of gene isoforms based on the PacBio IsoSeq method. Additionally, three individual follow-up analyses can be performed.

### Required: Isoform identification
An initial and always required step makes use of the [Iso-Seq](https://isoseq.how/) workflow to generate isoforms from PacBio long reads. In addition to the provided collapsing algorithm of PacBio, [TAMA](https://github.com/GenomeRIK/tama/wiki/Tama-Collapse) is also used as a second and more configurable algorithm, providing the "Wobble Walking" technique. Additionally, this step is able to perform demultiplexing of several SMRT-cell sequencing runs and combines duplicate samples present in different SMRT-cells. Note, that the clustering step in the Iso-Seq protocol is skipped to maximize the number of possible transcripts considered as individual isoform.

### Follow up 1: Sqanti
[SQANTI3](https://github.com/ConesaLab/SQANTI3) is a tool which compares the newly generated isoforms to an already existing annotation and classifies them in to predefined categories [Sqanti-categories](https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-isoform-classification:-categories-and-subcategories). Additionally it reports several output metrics ([Sqanti-Output](https://github.com/ConesaLab/SQANTI3/wiki/Understanding-the-output-of-SQANTI3-QC#classifcols)) in table format and reports several statistics in PDF format automatically.

### Follow up 2: PASA
The second follow-up analysis uses the tool [PASA](https://github.com/PASApipeline/PASApipeline/wiki), which uses the newly annotated isoforms to update an exsisting annotation.

### Follow up 3: Transdecoder
The tool transdecoder [Transdecoder](https://github.com/TransDecoder/TransDecoder/wiki) identifies potential coding regions within the transcript sequences. Note, that the ```--complete_orfs_only``` flag will be used per default.

# How to install
To run the pipeline, create a snakemake 7 environment: ```mamba create -n snakemake7 snakemake=7.22 python=3.11```\
Activate it using ```mamba activate snakemake7```\
All packages for the Isoform identification, PASA and Transdecoder will be installed automatically.\
For Sqanti,follow the official [installation instructions](https://github.com/ConesaLab/SQANTI3/wiki/Dependencies-and-installation).

Next, create the working directory you want to run the pipeline in and copy the following files to this directory:
- ```Snakefile```
- ```config_template.yml```
  
Then, copy the following folders with all files to the same directory:
- ```rules/```
- ```scripts/```
  
And last, specify the paths to required input files for the different steps of the pipeline.

# Config file
The main configurations of all inputs and the main analysis tools are defined via a config file (```config_template.yaml```). An example is provided in the ```config_template_example.yaml```.

### General setup
  - ```workdir```: All analyses will be done relative to the directory specified in ```workdir```. This directory also specifies the folder in which all required input files are stored.
  - ```pools```: defines the names of individual SMRT-cell runs in list format which is used in the workflow to uniquely identify input subread-files and demultiplexed outputs.
  - ```sample2pools```: is a dictionary of lists, mapping the sample IDs to the SMRT-cell run(s) in which they were pooled for sequencing.
  - ```name2bc```: is a list of dictionaries. Each dictionary represents a pool (in the order defined in ```pools```) and maps the sample ID to their repsective barcode used in the SMRT-cell run.

### General input files
  - ```ref_genomes```: path to reference genomes using the ```{sample}``` pattern specified in the general setup (required in all steps)
  - ```ref_annotations```: path to reference annotations. required patterns are ```{sample}``` (corresponding to the ones specified in the general setup as well as their reference genomes) and ```{tool}```, defining one or more techiques which were used to generate the respective reference annotation type. This allows to run follow up analyses with different reference analyses (used in follow up 1 and 2). Even if only one annotation type is present, a general pattern has to be chosen here.
  - ```barcodes```: barcode file, defining the sequences of 5' and 3' barcoded primers in fasta format. See the [PacBio CLI-workflow]https://isoseq.how/clustering/cli-workflow.html) for more details (required in the initial Isoform identification)
  - ```biosamples```: One file per ```{pool}```, mapping the primer IDs to the sample name in the SMRT-cell run. See the [Lima documentation](https://lima.how/faq/biosample.html) for more information (required in the initial Isoform identification). See ```input/meta/biosample_pool1_example.csv```.

### CCS
  - ```ccs_subreads```: path to the subreads. ```{pool}``` is required.
  - ```ccs_threads```: Threads to use during consensus calling
  - ```ccs_params```: paramters set (default none)

### Lima
  - ```lima_params```: parameters for the tool lima
  - ```lima_threads```: threads used for demultiplexing

### Refine
  - ```refine_params```: parameters for refine
  - ```refine_threads```: threads used for read refinement

### Isoform collaspsing (PacBio)
  - ```pbmm2```: PacBio minimap-wrapper parameters
  - ```pb_c```: PacBio collapsing parameters

### Isoform collaspsing (Tama)
  - ```minimap2```: minimap2 maping parameters (default: mimic pbmm2 parameters)
  - ```tama_c```: Tama collapsing parameters

### Isoform generation
```version```: If several collapsing parameters are tested (one at a time), a version can be specified (default "v1"). A subfolder will be created in the output directory for every version, in which all subsequent steps are stored.

### Sqanti
  - ```sqanti_path```: path to sqanti_qc.py
    
### PASA
  - ```pasa_set```: The parameter set used by PASA to integrate the isoforms into the reference annotations. Parameters are specified with ```pasa_params```
  - ```pasa_params```:TODO
  - ```config_align```: path/to/pasa.alignAssembly.Template.txt
  - ```config_compare```: path/to/pasa.annotationCompare.Template.txt

### Run
Specify the config file to use on top of the ```Snakefile``` and then start the pipeline with the following command:
```snakemake -j N --use-conda```, where N is the amount of jobs to run simultaneously
