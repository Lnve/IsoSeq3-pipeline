### SET PYTHONPATH!
#export PYTHONPATH=$PYTHONPATH:/ebio/abt6/lvaness/software/cDNA_Cupcake/sequence

rule sqanti_gff2gtf:
    input:
        "output/isoform_generation/collapse_{run}/{v}/{sample}.collapsed.{run}.{v}.gff"
    output:
        "output/isoform_generation/collapse_{run}/{v}/{sample}.collapsed.{run}.{v}.gtf"
    log:
        "log/sqanti/gff2gtf/{v}/{sample}.{run}.gff2gtf.{v}.log"
    conda:      
        "SQANTI3.env"
    benchmark:
        "stats/sqanti/gff2gtf/{v}/{sample}.{run}.gff2gtf.{v}.txt"
    shell:
        "gffread -T -o {output} {input}"

rule squanti:
    input:
        #annotation="../../input/annotations/{tool}/{sample}.{tool}-v2.3.gtf",
        annotation=config["ref_annotations"],
        genome=config['ref_genomes'],
        isoforms="output/isoform_generation/collapse_{run}/{v}/{sample}.collapsed.{run}.{v}.gtf"
    output:
        directory("output/sqanti/{run}/{v}/{sample}_{tool}_sqanti/")
    params:
        prefix="{sample}-{run}-vs-{tool}",
        settings="--report both --skipORF",
        path=config['sqanti_path']
    log:
        "log/sqanti/sqanti/{v}/{sample}.{run}.{tool}.sqanti.{v}.log"
    conda:
        "SQANTI3.env"
    benchmark:
        "stats/sqanti/sqanti/{v}/{sample}.{run}.{tool}.sqanti.{v}.txt"
    shell:
        "export PYTHONPATH=\$PYTHONPATH:/ebio/abt6/lvaness/software/cDNA_Cupcake/sequence && mkdir -p {output} && python3 {params.path} --report pdf --skipORF -o {params.prefix} --d {output} {input.isoforms} {input.annotation} {input.genome} > {log}"
