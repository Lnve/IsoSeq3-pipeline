##############################
###### Read Preparation ######
##############################

#-------------#
##### CCS #####
#-------------#

rule ccs:
    """
    Create HiFi reads from SequelII IsoSeq subreads
    """
    input:
        config["ccs_subreads"]
    output:
        "output/read_preparation/ccs/{pool}/{pool}.hifi.bam"
    log:
        "log/read_preparation/ccs/{pool}.ccs.log"
    params:
        outdir="output/read_preparation/ccs/{pool}",
        prefix="output/read_preparation/ccs/{pool}/{pool}",
        flags=config['ccs_params']
    threads: config['ccs_threads']
    conda:
        "envs/isoseq.yml"
    benchmark:
        "stats/read_preparation/ccs/{pool}.ccs.txt"
    shell:
        "mkdir -p {params.outdir} && ccs --report-file {params.prefix}.report.txt --metrics-json {params.prefix}.metrics.json.gz --hifi-summary-json {params.prefix}.hifi-summary.json {params.flags} --num-threads {threads} --log-level INFO --log-file {log} {input} {output}"

#---------------#
##### Lima ######
#---------------#

for idx, pool in enumerate(config['pools']):
    n2bc = config['name2bc'][idx]
    bcs = list(n2bc.values())
    rule lima:
        """
        Demultiplex the pools and remove primers + barcodes (TO DO: do the
        barcodes present also contain primers? --> Barcodes + universal primers
        from the downloaded isoseq file)
        """
        name: f"lima_{pool}"
        input:
            hifi=f"output/read_preparation/ccs/{pool}/{pool}.hifi.bam",
            biosamples=config['biosamples'].format(pool=pool),
            barcodes=config['barcodes']
        output: 
            expand(f"output/read_preparation/lima/{pool}/{pool}.demux.hifi.{{bc}}.bam", bc=bcs )
        log:
            f"log/read_preparation/lima/{pool}.lima.log"
        params:
            pattern=f"output/read_preparation/lima/{pool}/{pool}.demux.hifi.bam",
            flags=config['lima_params']
        threads: config['lima_threads']
        conda:
            "envs/isoseq.yml"
        benchmark:
            f"stats/read_preparation/lima/{pool}.lima.txt"
        shell:
            "lima"
            " {params.flags}"
            " --num-threads {threads}"
            " --log-level INFO"
            " --log-file {log}"
            " --store-unbarcoded"
            " --split-bam"
            " --biosample-csv {input.biosamples}"
            " {input.hifi}"
            " {input.barcodes}"
            " {params.pattern}"

    #rule renameit:
    #    name: f"renameit_{pool}_{sample}"
    #    input:
    #        lambda wc: f"output/read_preparation/lima/{pool}/{pool}.demux.hifi.{n2bc[wc.sample]}.bam"
    #        #f"output/read_preparation/lima/{pool}/{pool}.demux.hifi.{b}.bam"
    #    output:
    #        f"output/read_preparation/lima/{pool}/{pool}.demux.hifi.{{sample}}.renamed.bam"
    #    shell:
    #        "ln -sf $(readlink -f {input}) {output}"

    for b in bcs:
        rule renameit:
            #name: f"renameit_{pool}_{sample}"
            name: f"renameit_{pool}"
            input:
                f"output/read_preparation/lima/{pool}/{pool}.demux.hifi.{b}.bam"
            output:
                f"output/read_preparation/lima/{pool}/{pool}.demux.hifi.{{sample}}.renamed.bam"
            shell:
                "ln -sf $(readlink -f {input}) {output}"


#----------------#
##### Refine #####
#----------------#

rule refine:
    """
    Remove polyAtails and filter for concatamers (primers required again)
    """
    input:
        reads="output/read_preparation/lima/{pool}/{pool}.demux.hifi.{sample}.renamed.bam",
        primers=config['barcodes']
    output:
        "output/read_preparation/refine/{pool}/{pool}.{sample}.ncfl.bam"
    log:
        "log/read_preparation/refine/{pool}.{sample}.refine.log"
    params:
        dir=directory("output/read_preparation/refine/{pool}"),
        flags=config['refine_params']
    threads: config['refine_threads']
    conda:
        "envs/isoseq.yml"
    benchmark:
        "stats/read_preparation/refine/{pool}.{sample}.refine.txt"
    shell:
        "mkdir -p {params.dir} && isoseq3 refine -j {threads} {params.flags} --log-level INFO --log-file {log} {input.reads} {input.primers} {output}"

#---------------------------#
##### Rename read group #####
#---------------------------#

rule rename_rg:
    input:
        "output/read_preparation/refine/{pool}/{pool}.{sample}.ncfl.bam"
    output:
        temp("output/read_preparation/rename_rg/{pool}/{pool}.{sample}.ncfl.bam")
    log:
        "log/read_preparation/rename_rg/{pool}.{sample}.rename_rg.log"
    params:
        dir=directory("output/read_preparation/rename_rg/{pool}"),
        pattern="{pool}_{sample}"
    threads: 4
    conda:
        "envs/isoseq.yml"
    benchmark:
        "stats/combine/{pool}_{sample}_rename_rg_benchmark.txt"
    shell:
        "mkdir -p {params.dir} && (python3 scripts/rg.py -t {threads} -i {input} -o {output} -d {params.pattern}) > {log}"

#------------------------------#
##### Concatenate all runs #####
#------------------------------#

from collections import defaultdict

rule concat_samples:
    input:
        lambda wc: expand("output/read_preparation/rename_rg/{pool}/{pool}.{sample}.ncfl.bam",
                    pool=config["sample2pools"][wc.sample],
                    sample=wc.sample)
    output:
        "output/read_preparation/concat_samples/{sample}.ncfl.bam"
    log:
        "log/read_preparation/concat_samples/{sample}.concat_samples.log"
    conda:
        "envs/isoseq.yml"
    benchmark:
        "stats/read_preparation/concat_samples/{sample}.concat_samples.txt"
    shell:
        "mkdir -p output/read_preparation/concat_samples && (samtools merge {output} {input}) > {log}"
