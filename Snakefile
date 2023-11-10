configfile: "config_template.yaml"

workdir : config["workdir"]

# general parameters:
v = config["version"]

ruleorder:  sqanti_gff2gtf > transcdecoder_gff2gtf


refine_out = []
rename_rg_out = []
concat_out = []
cluster_out = []


for idx, pool in enumerate(config['pools']):
    for bc, sample in zip(config['name2bc'][idx].values(), config['name2bc'][idx].keys()):
        refine_out.append(f"{pool}/{pool}.{sample}.ncfl.bam")
        rename_rg_out.append(f"{pool}/{pool}.{sample}.ncfl.bam")
        concat_out.append(f"{sample}.ncfl.bam")
        cluster_out.append(f"{sample}.transcripts.bam")

# read preparation
include: "rules/read_preparation.smk"

# isoform collapsing
include: "rules/isoform_generation.smk"

# sqanti
include: "rules/sqanti.smk"

# pasa
include: "rules/pasa.smk"

# transdecoder
include: "rules/transdecoder.smk"

rule all:
    input:
        # read preparation
        #expand("output/read_preparation/ccs/{pool}/{pool}.hifi.bam", pool=config["pools"]),
        #[f"output/read_preparation/lima/{pool}/{pool}.demux.hifi.{sample}.renamed.bam"
        #    for sample, pools in config["sample2pools"].items()
        #    for pool in pools],
        #expand("output/read_preparation/refine/{pattern}", pattern=refine_out),
        #expand("output/read_preparation/rename_rg/{pattern}", pattern=rename_rg_out),
        expand("output/read_preparation/concat_samples/{pattern}", pattern=concat_out),
        
        ### collapsing.smk
        expand("output/isoform_generation/collapse_pb/{v}/{sample}.collapsed.pb.{v}.gff", sample=config["samples2map"], v=v),
        expand("output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}.gff3", sample=config["samples2map"], v=v),
        expand("output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}_read_support.txt", sample=config["samples2map"], v=v),
        
        ### sqanti.smk
        #expand("output/sqanti/{run}/{v}/{sample}_{tool}_sqanti/", sample=config["samples_sqanti"], v=v, run=["pb", "tama"], tool=["liftoff","augustus"]),
        
        ### pasa.smk
        expand("output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/compare1.done", sample=config["samples_pasa"], v=v, run=["tama"], params=config['pasa_settings'], tool=["liftoff", "augustus"]),
        
        ### Transdecoder
        #expand("output/transdecoder/{run}/{v}/{sample}_transdecoder/{sample}_transcripts.fasta.transdecoder.genome.sorted.gff3", sample=config["samples_transdecoder"], v=v, run=["tama", "pb"]),
