# IsoSeq3-pipeline
The pipeline provides an automated Snakemake workflow for the identification and annotation of gene isoforms based on the PacBio IsoSeq method. Additionally, three individual follow-up analyses can be performed.

### Required: Isoform identification
An initial and always required step makes use of the [Iso-Seq](https://isoseq.how/) workflow to generate isoforms from PacBio long reads. In addition to the provided collapsing algorithm of PacBio, [TAMA](https://github.com/GenomeRIK/tama/wiki/Tama-Collapse) is also used as a second and more configurable algorithm, providing the "Wobble Walking" technique. Additionally, this step is able to perform demultiplexing of several SMRT-cell sequencing runs and combines duplicate samples present in different SMRT-cells. Note, that the clustering step in the Iso-Seq protocol is skipped to maximize the number of possible transcripts considered as individual isoform.

### Follow up 1: Sqanti
[SQANTI3](https://github.com/ConesaLab/SQANTI3) is a tool which compares the newly generated isoforms to an already existing annotation and classifies them in to predefined categories [Sqanti-categories](https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-isoform-classification:-categories-and-subcategories). Additionally it reports several output metrics ([Sqanti-Output](https://github.com/ConesaLab/SQANTI3/wiki/Understanding-the-output-of-SQANTI3-QC#classifcols)) in table format and reports several statistics in PDF format automatically.

### Follow up 2: PASA
The second follow-up analysis uses the tool [PASA](https://github.com/PASApipeline/PASApipeline/wiki), which uses the newly annotated isoforms to update an exsisting annotation.

### Follow up 3: Transdecoder
The tool transdecoder [Transdecoder](https://github.com/TransDecoder/TransDecoder/wiki) identifies potential coding regions within the transcript sequences. Note, that the '''--complete_orfs_only''' flag will be used per default.

### How to install
TODO

### Config file
The main configurations of all inputs and the main analysis tools are defined via a config file ('config_template.yaml'). 

### Run 
