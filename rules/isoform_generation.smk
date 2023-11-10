#################################
########## Preparation ##########
#################################

#rule bam2fastq:
#    input:
#        "output/read_preparation/concat_samples/{sample}.ncfl.bam"
#    output:
#        "output/read_preparation/concat_samples/{sample}.ncfl.fastq"
#    log:
#        "log/isoform_generation/bam2fastq/{sample}.bam2fastq.log"
#    conda:
#        "isoseq"
#    benchmark:
#        "stats/isoform_generation/bam2fastq/{sample}.bam2fastq.txt"
#    shell:
#        "bedtools bamtofastq -i {input} -fq {output} > {log}" #doubles reads
#
#rule rm_dups:
#    input:
#        "output/read_preparation/concat_samples/{sample}.ncfl.fastq"
#    output:
#        "output/read_preparation/concat_samples/{sample}.ncfl.no_dups.fastq"
#    log:
#        "log/isoform_generation/rm_dups/{sample}.rm_dups.log"
#    params:
#        stats="output/concat_samples/{sample}.ncfl.dups.stats"
#    priority: 50
#    threads: 4
#    conda:
#        "isoseq"
#    benchmark:
#        "stats/isoform_generation/rm_dups/{sample}.rm_dups.txt"
#    shell:
#        "(cat {input} | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n'  | seqkit rmdup -n -j {threads} -D {params.stats} -o {output}) > tee {log}"

rule bam2fastq:
    input:
        "output/read_preparation/concat_samples/{sample}.ncfl.bam"
    output:
        "output/isoform_generation/bam2fastq/{sample}.ncfl.fastq"
    log:
        "log/isoform_generation/bam2fastq/{sample}.bam2fastq.log"
    conda:
        "envs/isoseq.yml"
    benchmark:
        "stats/isoform_generation/bam2fastq/{sample}.bam2fastq.txt"
    shell:
        "(samtools bam2fq {input} > {output}) > {log}"


##########################
########### PB ###########
##########################

########## Mapping ##########
rule pbmm2:
    input:
        ref=config['ref_genomes'],
        transcripts="output/isoform_generation/bam2fastq/{sample}.ncfl.fastq"
    output:
        "output/isoform_generation/pbmm2/{sample}.mapped.bam"
    log:
        "log/isoform_generation/pbmm2/{sample}.pbmm2.log"
    params:
        pbmm2=config['pbmm2']
    threads: 10
    conda:
        "envs/isoseq.yml"
    benchmark:
        "stats/isoform_generation/pbmm2/{sample}.pbmm2.txt"
    shell:
        "mkdir -p output/isoform_generation/pbmm2 && pbmm2 align {params.pbmm2} --log-level INFO --log-file {log} {input.transcripts} {input.ref} {output}"

########## Collapsing ##########

rule collapse_pb:
    input:
        "output/isoform_generation/pbmm2/{sample}.mapped.bam"
    output:
        "output/isoform_generation/collapse_pb/{v}/{sample}.collapsed.pb.{v}.gff"
    log:
        "log/isoform_generation/collapse_pb/{v}/{sample}.collapse_pb.{v}.log"
    threads: 4
    params:
        v=config["pb_c"]
    conda:
        "envs/isoseq.yml"
    benchmark:
        "stats/isoform_generation/collapse_pb/{v}/{sample}.collapse_pb.{v}.txt"
    shell:
        "mkdir -p output/isoform_generation/collapse_pb/{v} && isoseq3 collapse {params.v} --log-level INFO --log-file {log} -j {threads} {input} {output}"


##########################
########## TAMA ##########
##########################

########## Mapping ##########

rule minimap2:
    """
    Map the reads against the genome. Don't use pbmm2, as I don't know if -ax splice is enabled here, which is required for StringTiea, BUT use the equivalent parameters as ISOSEQ preset from pbmm2 (https://github.com/lh3/minimap2/issues/769).
    """
    input:
        ref=config['ref_genomes'],
        transcripts="output/isoform_generation/bam2fastq/{sample}.ncfl.fastq"
    output:
        "output/isoform_generation/minimap2/{sample}.mapped.bam"
    log:
        "log/isoform_generation/minimap2/{sample}.minimap2.log"
    params:
        minimap2=config['minimap2']
    threads: 10
    conda:
        "envs/isoseq.yml"
    benchmark:
        "stats/isoform_generation/minimap2/{sample}.minimap2.txt"
    shell:
        "mkdir -p output/isoform_generation/minimap2 && (minimap2 -t {threads} {params.minimap2} {input.ref} {input.transcripts} | samtools sort -l 9 -o {output}) > {log}"

rule minimap2_index:
    input:
        "output/isoform_generation/minimap2/{sample}.mapped.bam"
    output:
        "output/isoform_generation/minimap2/{sample}.mapped.bam.bai"
    log:
        "log/isoform_generation/minimap2_index/{sample}.minimap2_index.log"
    threads: 10
    conda:
        "envs/isoseq.yml"
    benchmark:
        "stats/isoform_generation/minimap2_index/{sample}.minimap2_index.txt"
    shell:
        "samtools index -o {output} {input} > {log}"

########## Collapsing ##########

rule collapse_tama:
    input:
        mapping="output/isoform_generation/minimap2/{sample}.mapped.bam",
        fasta=config['ref_genomes']
    output:
        "output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}_collapsed.bed"
    log:
        "log/isoform_generation/collapse_tama/{v}/{sample}.collapse_tama.{v}.log"
    params:
        python_path="~/anaconda3/envs/py2.7/bin/python2.7",
        tama_path="~/software/tama/tama_collapse.py",
        prefix="output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}",
        v=config["tama_c"]
    conda:
        "envs/tama.yml"
    benchmark:
        "stats/isoform_generation/collapse_tama/{v}/{sample}.collapse_tama.{v}.txt"
    shell:
        # "mkdir -p output/{v}/collapse_tama && ({params.python_path} {params.tama_path} -s {input.mapping} -b BAM -f {input.fasta} -p {params.prefix} {params.v} -log log_off) > {log}"
        # "mkdir -p output/isoform_generation/collapse_tama/{v} && (python2.7 {params.tama_path} -s {input.mapping} -b BAM -f {input.fasta} -p {params.prefix} {params.v} -log log_off) > {log}"
        "mkdir -p output/isoform_generation/collapse_tama/{v} && (tama_collapse.py -s {input.mapping} -b BAM -f {input.fasta} -p {params.prefix} {params.v} -log log_off) > {log}"

rule bed2gtf:
    input:
        "output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}_collapsed.bed"
    output:
        "output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}.gtf"
    log:
        "log/isoform_generation/bed2gtf/{v}/{sample}.bed2gtf.{v}.log"
    benchmark:
        "stats/isoform_generation/bed2gtf/{v}/{sample}.bed2gtf.{v}.txt"
    shell:
        "(perl scripts/bed12ToGTF.1.pl < {input} > {output}) > {log}"

rule gtf2gff3:
    input:
        "output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}.gtf"
    output:
        "output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}.gff3"
    log:
        "log/isoform_generation/gtf2gff3/{v}/{sample}.gtf2gff3.{v}.log"
    conda:
        "envs/isoseq.yml"
    benchmark:
        "stats/isoform_generation/gtf2gff3/{v}/{sample}.gtf2gff3.{v}.txt"
    shell:
        "gffread -E --keep-genes {input} -o {output} > {log}"

######### Read Support #########

rule tama_file_list:
    input:
        "output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}.gff3"
    output:
        "output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}.reads.txt"
    log:
        "log/isoform_generation/tama_file_list/{v}/{sample}.tama_file_list.{v}.log"
    conda:
        "envs/isoseq.yml"
    benchmark:
        "stats/isoform_generation/tama_file_list/{v}/{sample}.tama_file_list.{v}.txt"
    shell:
        """
        touch {input}
        (echo -e 'flnc\toutput/isoform_generation/collapse_tama/{wildcards.v}/{wildcards.sample}.collapsed.tama.{wildcards.v}_trans_read.bed\ttrans_read' > {output}) > {log}
        """

rule tama_read_support:
    input:
        "output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}.reads.txt"
    output:
        "output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}_read_support.txt"
    log:
        "log/isoform_generation/tama_read_support/{v}/{sample}.tama_read_support.{v}.log"
    params:
        tama_path="~/software/tama/tama_go/read_support/tama_read_support_levels.py",
        prefix="output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}"
    conda:
        "envs/tama.yml"
    benchmark:
        "stats/isoform_generation/tama_read_support/{v}/{sample}.tama_read_support.{v}.txt"
    shell:
        # "python {params.tama_path} -f {input} -o {params.prefix} -m no_merge > {log}"
        "tama_read_support_levels.py -f {input} -o {params.prefix} -m no_merge > {log}"
