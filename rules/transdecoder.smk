rule transcdecoder_gff2gtf:
    input:
        "output/isoform_generation/collapse_{run}/{v}/{sample}.collapsed.{run}.{v}.gff"
    output:
        "output/isoform_generation/collapse_{run}/{v}/{sample}.collapsed.{run}.{v}.gtf"
    log:
        "log/transdecoder/gff2gtf/{v}/{sample}.{run}.gff2gtf.{v}.log"
    conda:      
        "envs/isoseq.yml"
    benchmark:
        "stats/sqanti/gff2gtf/{v}/{sample}.{run}.gff2gtf.{v}.txt"
    shell:
        "gffread -T -o {output} {input}"

rule add_ID_tama:
    input:
        "output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}.gtf"
    output:
        "output/transdecoder/tama/{v}/{sample}.collapsed.tama.{v}.addID.gtf"
    log:
        "log/transdecoder/add_ID_tama/{v}/{sample}.add_ID_tama.{v}.log"
    conda:
        "envs/pasa.yml"
    benchmark:
        "stats/transdecoder/add_ID_tama/{v}/{sample}.add_ID_tama.{v}.txt"
    shell:
        "python3 scripts/transdecoder_conversion.py -i {input} -o {output} > {log}"

rule gtf2fasta_tama:
    input:
        gtf="output/transdecoder/tama/{v}/{sample}.collapsed.tama.{v}.addID.gtf",
        fasta=config['ref_genomes']
    output:
        fasta="output/transdecoder/tama/{v}/{sample}_transdecoder/{sample}_transcripts.fasta",
        check=temp("output/transdecoder/tama/{v}/{sample}_transdecoder/fasta.done")
    log:
        "log/transdecoder/gtf2fasta_tama/{v}/{sample}.gtf2fasta_tama.{v}.log"
    params:
        out=directory("output/transdecoder/tama/{v}/{sample}_transdecoder/")
    conda:
        "envs/pasa.yml"
    benchmark:
        "stats/transdecoder/gtf2fasta_tama/{v}/{sample}.gtf2fasta_tama.{v}.txt"
    shell:
        "mkdir -p {params.out} && (gtf_genome_to_cdna_fasta.pl {input.gtf} {input.fasta} > {output.fasta}) > {log} && touch {output.check}"

rule gtf2fasta_pb:
    input:
        gtf="output/isoform_generation/collapse_pb/{v}/{sample}.collapsed.pb.{v}.gtf",
        fasta=config['ref_genomes']
    output:
        fasta="output/transdecoder/pb/{v}/{sample}_transdecoder/{sample}_transcripts.fasta",
        check=temp("output/transdecoder/pb/{v}/{sample}_transdecoder/fasta.done")
    log:
        "log/transdecoder/gtf2fasta_pb/{v}/{sample}_gtf2fasta_pb.{v}.log"
    params:
        out=directory("output/transdecoder/pb/{v}/{sample}_transdecoder/")
    conda:
        "envs/pasa.yml"
    benchmark:
        "stats/transdecoder/gtf2fasta_pb/{v}/{sample}_gtf2fasta_pb.{v}.txt"
    shell:
        "mkdir -p {params.out} && (gtf_genome_to_cdna_fasta.pl {input.gtf} {input.fasta} > {output.fasta}) > {log} && touch {output.check}"

rule gtf2gff_tama:
    input:
        "output/transdecoder/tama/{v}/{sample}.collapsed.tama.{v}.addID.gtf"
    output:
        "output/transdecoder/tama/{v}/{sample}_transdecoder/{sample}_transcripts.gff3"
    log:
        "log/transdecoder/gtf2gff_tama/{v}/{sample}.gtf2gff_tama.{v}.log"
    params:
        path="output/{v}/collapse_tama/{sample}_transdecoder/"
    conda:
        "envs/pasa.yml"
    benchmark:
        "stats/transdecoder/gtf2gff_tama/{v}/{sample}.gtf2gff_tama.{v}.txt"
    shell:
        "(gtf_to_alignment_gff3.pl {input} > {output}) > {log}"

rule gtf2gff_pb:
    input:
        "output/isoform_generation/collapse_pb/{v}/{sample}.collapsed.pb.{v}.gtf"
    output:
        "output/transdecoder/pb/{v}/{sample}_transdecoder/{sample}_transcripts.gff3"
    log:
        "log/transdecoder/gtf2gff_pb/{v}/{sample}.gtf2gff_pb.{v}.log"
    params:
        path="output/{v}/collapse_ob/{sample}_transdecoder/"
    conda:
        "envs/pasa.yml"
    benchmark:
        "stats/transdecoder/gtf2gff_pb/{v}/{sample}.gtf2gff_pb.{v}.txt"
    shell:
        "(gtf_to_alignment_gff3.pl {input} > {output}) > {log}"

rule longORF:
    input:
        check="output/transdecoder/{run}/{v}/{sample}_transdecoder/fasta.done"
    output:
        temp("output/transdecoder/{run}/{v}/{sample}_transdecoder/longORF.done")
    log:
        "log/transdecoder/longORF/{v}/{sample}.{run}.longORF.{v}.log"
    params:
        path="output/transdecoder/{run}/{v}/{sample}_transdecoder/",
        fasta="{sample}_transcripts.fasta"
    conda:
        "envs/pasa.yml"
    benchmark:
        "stats/transdecoder/longORF/{v}/{sample}.{run}.longORF.{v}.txt"
    shell:
        """
        touch {input.check}
        cd {params.path}
        TransDecoder.LongOrfs -t {params.fasta} > ../../../../../{log}
        cd ../../../../../
        touch {output}
        """

rule predictORF:
    input:
        check="output/transdecoder/{run}/{v}/{sample}_transdecoder/longORF.done"
    output:
        temp("output/transdecoder/{run}/{v}/{sample}_transdecoder/Predict.done")
    log:
        "log/transdecoder/predictORF/{v}/{sample}.{run}.predictORF.{v}.log"
    params:
        path="output/transdecoder/{run}/{v}/{sample}_transdecoder/",
        fasta="{sample}_transcripts.fasta"
    conda:
        "envs/pasa.yml"
    benchmark:
        "stats/transdecoder/predictORF/{v}/{sample}.{run}.predictORF.{v}.txt"
    shell:
        """
        touch {input.check}
        cd {params.path}
        TransDecoder.Predict -t {params.fasta} > ../../../../../{log}
        cd ../../../../../
        touch {output}
        """

rule genome_anno:
    input:
        check="output/transdecoder/{run}/{v}/{sample}_transdecoder/Predict.done",
        gff3="output/transdecoder/{run}/{v}/{sample}_transdecoder/{sample}_transcripts.gff3",
        fasta="output/transdecoder/{run}/{v}/{sample}_transdecoder/{sample}_transcripts.fasta"
    output:
        "output/transdecoder/{run}/{v}/{sample}_transdecoder/{sample}_transcripts.fasta.transdecoder.genome.gff3"
    log:
        "log/transdecoder/genome_anno/{v}/{sample}.{run}.genome_anno.{v}.log"
    params:
        trans_gff="output//transdecoder/{run}/{v}/{sample}_transdecoder/{sample}_transcripts.fasta.transdecoder.gff3"
    conda:
        "envs/pasa.yml"
    benchmark:
        "stats/transdecoder/genome_anno/{v}/{sample}.{run}.genome_anno.{v}.txt"
    shell:
        "(touch {input.check} && cdna_alignment_orf_to_genome_orf.pl {params.trans_gff} {input.gff3} {input.fasta} > {output}) > {log}"

rule transdecoder_clean:
    input:
        "output/transdecoder/{run}/{v}/{sample}_transdecoder/{sample}_transcripts.fasta.transdecoder.genome.gff3"
    output:
        "output/transdecoder/{run}/{v}/{sample}_transdecoder/{sample}_transcripts.fasta.transdecoder.genome.cleaned.gff3"
    log:
        "log/transdecoder/transdecoder_clean/{v}/{sample}.{run}.transdecoder_clean.{v}.log"
    conda:
        "envs/pasa.yml"
    benchmark:
        "stats/transdecoder/transdecoder_clean/{v}/{sample}.{run}.transdecoder_clean.{v}.txt"
    shell:
        "(cat {input} | grep {wildcards.sample} > {output}) > {log}"

rule transdecoder_sort:
    input:
        "output/transdecoder/{run}/{v}/{sample}_transdecoder/{sample}_transcripts.fasta.transdecoder.genome.cleaned.gff3"
    output:
        "output/transdecoder/{run}/{v}/{sample}_transdecoder/{sample}_transcripts.fasta.transdecoder.genome.sorted.gff3"
        # "/ebio/abt6_projects7/diffLines_20/data/Luisa_web_apollo/20230202_new_files_for_webapollo_v2/20230202_new_files_for_webapollo/transdecoder2/{sample}.{run}{v}_transcripts.fasta.transdecoder.genome.sorted.gff3"
    log:
        "log/transdecoder/transdecoder_sort/{v}/{sample}.{run}.transdecoder_sort.{v}.log"
    conda:
        "envs/pasa.yml"
    benchmark:
        "stats/transdecoder/transdecoder_sort/{v}/{sample}.{run}.transdecoder_sort.{v}.txt"
    shell:
        "agat_convert_sp_gxf2gxf.pl -g {input} -o {output} > {log}"
