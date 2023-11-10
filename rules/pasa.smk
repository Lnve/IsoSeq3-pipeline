rule generate_iso_fasta_tama:
    input:
        genome=config['ref_genomes'],
        gff="output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}.gff3"
    output:
        fasta="output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}.fasta"
    log:
        "log/pasa/generate_iso_fasta_tama/{v}/{sample}.tama.generate_iso_fasta_tama.{v}.log"
    priority: 100
    conda:
        "envs/isoseq.yml"
    benchmark:
        "stats/pasa/generate_iso_fasta_tama/{v}/{sample}.tama.generate_iso_fasta_tama.txt"
    shell:
        "gffread -w {output.fasta} -g {input.genome} {input.gff} > {log}"

rule generate_iso_fasta_pb:
    input:
        genome="config['ref_genomes']",
        gff="output/isoform_generation/collapse_pb/{v}/{sample}.collapsed.pb.{v}.gff"
    output:
        fasta="output/isoform_generation/collapse_pb/{v}/{sample}.collapsed.pb.{v}.fasta"
    log:
        "log/pasa/generate_iso_fasta_pb/{v}/{sample}.pb.generate_iso_fasta_pb.{v}.log"
    priority: 100
    conda:
        "envs/isoseq.yml"     
    benchmark:
        "stats/pasa/generate_iso_fasta_pb/{v}/{sample}.pb.generate_iso_fasta_pb.{v}.txt"
    shell:
        "gffread -w {output.fasta} -g {input.genome} {input.gff} > {log}"

rule seqclean:
    input:
        "output/isoform_generation/collapse_{run}/{v}/{sample}.collapsed.{run}.{v}.fasta"
    output:
        "output/pasa/{run}/{v}/{sample}_pasa_seqclean.done"
    log:
        config['workdir'] + "/log/pasa/seqclean/{v}/{sample}.{run}.seqclean.{v}.log"
    params:
        dir=directory("output/isoform_generation/collapse_{run}/{v}/"),
        input="{sample}.collapsed.{run}.{v}.fasta"
    threads: 10
    conda:             
        "envs/pasa.yml"
    priority: 90
    benchmark:
        "stats/pasa/seqclean/{v}/{sample}.{run}.seqclean.{v}.txt"
    shell:
        """
        touch {input}
        cd {params.dir}
        pwd
        # ~/software/PASApipeline/bin/seqclean {params.input} -c {threads} > {log}
        $PASAHOME/bin/seqclean {params.input} -c {threads} > {log}
        cd ../../../../
        touch {output}
        """

rule config_align:
    input:
        config['config_align']
    output:
        cfg="output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/alignAssembly.config",
        check="output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/alignAssembly.done"
    log:
        "log/pasa/config_align/{v}/{sample}.{run}.{tool}.pasa_{params}.config_align.{v}.log"
    params:
        path=config['workdir'] + "/output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/{sample}_mydb_pasa",
        dir=directory("output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/")
    priority: 80
    conda:      
        "envs/pasa.yml"
    benchmark:
        "stats/pasa/config_align/{v}/{sample}.{run}.{tool}.pasa_{params}.config_align.{v}.txt"
    shell:
        """
        mkdir -p {params.dir}
        cp {input} {output.cfg}
        sed -i \"s@<__DATABASE__>@{params.path}@\" {output.cfg}
        touch {output.check}
        """

rule set_params:
    input:
        config['config_compare']
    output:
        "output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/annotCompare.config"
    log:
        "log/pasa/set_params/{v}/{sample}.{run}.{tool}.pasa_{params}.set_params.{v}.log"
    params:
        ov="30",
        rp="40",
        pp="70",
        fl="30",
        nfl="30",
        pa="30",
        pog="30",
        ue="2",
        dir=directory("output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/"),
    priority: 75
    conda:      
        "envs/pasa.yml"
    benchmark:
        "stats/pasa/set_params/{v}/{sample}.{run}.{tool}.pasa_{params}.set_params.{v}.txt"
    shell: 
        """
        mkdir -p {params.dir}
        python3 scripts/pasa_config.py -i {input} -o {output} -ov {params.ov} -rp {params.rp} -pp {params.pp} -fl {params.fl} -nfl {params.nfl} -pa {params.pa} -pog {params.pog} -ue {params.ue}
        """

rule config_anno:
    input:
        "output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/annotCompare.config"
    output:
        check="output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/annotCompare.done"
    log:
        "log/pasa/config_anno/{v}/{sample}.{run}.{tool}.pasa_{params}.config_anno.{v}.log"
    params:
        path=config['workdir'] + "/output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/{sample}_mydb_pasa",
        file="output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/annotCompare.config"
    priority: 80
    conda:      
        "envs/pasa.yml"
    benchmark:
        "stats/pasa/config_anno/{v}/{sample}.{run}.{tool}.pasa_{params}.config_anno.{v}.txt"
    shell:
        """
        sed -i \"s@<__DATABASE__>@{params.path}@\" {params.file}
        touch {output.check}
        """

rule pasa_alignment_tama:
    """
    run this prior to pasa_alignment_pb, since parallel indexing causes trouble
    """
    input:
        seqclean="output/pasa/tama/{v}/{sample}_pasa_seqclean.done",
        cfg="output/pasa/tama/{v}/{sample}_{tool}_pasa_{params}/alignAssembly.done"
    output:
        check="output/pasa/tama/{v}/{sample}_{tool}_pasa_{params}/alignment.done",
        idx_check="output/pasa/tama/{v}/{sample}_{tool}_pasa_{params}/idx_check.done"
    log:
        config['workdir'] + "log/pasa/pasa_alignment_tama/{v}/{sample}.tama.{tool}.pasa_{params}.pasa_alignment_tama.{v}.log"
    params:
        dir="output/pasa/tama/{v}/{sample}_{tool}_pasa_{params}/",
        config="alignAssembly.config",
        genome=config['workdir'] + "/" + config['ref_genomes'],
        fasta=config['workdir'] + "/output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}.fasta",
        fasta_clean=config['workdir'] + "/output/isoform_generation/collapse_tama/{v}/{sample}.collapsed.tama.{v}.fasta.clean",
        aligners="blat,gmap,minimap2"
    threads: 30
    priority: 70
    conda:      
        "envs/pasa.yml"
    benchmark:
        "stats/pasa/pasa_alignment_tama/{v}/{sample}.tama.{tool}.pasa_{params}.pasa_alignment_tama.{v}.txt"
    shell:
        """
        touch {input.seqclean}
        touch {input.cfg}
        cd {params.dir}
        # ~/software/PASApipeline/Launch_PASA_pipeline.pl -c {params.config} -C -R -g {params.genome} -t {params.fasta_clean} -T -u {params.fasta} --ALIGNERS {params.aligners} --CPU {threads} |& tee {log}
        ($PASAHOME/Launch_PASA_pipeline.pl -c {params.config} -C -R -g {params.genome} -t {params.fasta_clean} -T -u {params.fasta} --ALIGNERS {params.aligners} --CPU {threads}) > {log}
        cd ../../../../../
        touch {output.check}
        touch {output.idx_check}
        """

rule pasa_alignment_pb:
    input:
        seqclean="output/pasa/pb/{v}/{sample}_pasa_seqclean.done",
        cfg="output/pasa/pb/{v}/{sample}_{tool}_pasa_{params}/alignAssembly.done",
        idx_check="output/pasa/tama/{v}/{sample}_{tool}_pasa_{params}/idx_check.done"
    output:
        "output/pasa/pb/{v}/{sample}_{tool}_pasa_{params}/alignment.done"
    log:
        config['workdir'] + "log/pasa/pasa_alignment_pb/{v}/{sample}.pb.{tool}.pasa_{params}.pasa_alignment_pb.{v}.log"
    params:
        dir="output/pasa/pb/{v}/{sample}_{tool}_pasa_{params}/",
        config="alignAssembly.config",
        genome=config['workdir'] + "/" + config['ref_genomes'],
        fasta=config['workdir'] + "/output/isoform_generation/collapse_pb/{v}/{sample}.collapsed.pb.{v}.fasta",
        fasta_clean=config['workdir'] + "/output/isoform_generation/collapse_pb/{v}/{sample}.collapsed.pb.{v}.fasta.clean",
        aligners="blat,gmap,minimap2"
    threads: 30
    priority: 65
    conda:      
        "envs/pasa.yml"
    benchmark:
        "stats/pasa/pasa_alignment_pb/{v}/{sample}.pb.{tool}.pasa_{params}.pasa_alignment_pb.{v}.txt"
    shell:
        """
        touch {input.seqclean}
        touch {input.cfg}
        touch {input.idx_check}
        cd {params.dir}
        ($PASAHOME/Launch_PASA_pipeline.pl -c {params.config} -C -R -g {params.genome} -t {params.fasta_clean} -T -u {params.fasta} --ALIGNERS {params.aligners} --CPU {threads}) > {log}
        cd ../../../../../
        touch {output}
        """

rule pasa_load_anno1:
    input:
        "output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/alignment.done"
    output:
        "output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/load1.done"
    log:
        config['workdir'] + "log/pasa/pasa_load_anno1/{v}/{sample}.{run}.{tool}.pasa_{params}.pasa_load_anno1.{v}.log"
    params:
        dir="output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/",
        config="alignAssembly.config",
        genome=config['workdir'] + "/" + config['ref_genomes'],
        anno2update=config['workdir'] + "/" + config["ref_annotations"]
    priority: 60
    conda:      
        "envs/pasa.yml"
    benchmark:
        "stats/pasa/pasa_load_anno1/{v}/{sample}.{run}.{tool}.pasa_{params}.pasa_load_anno1.{v}.txt"
    shell:
        """
        touch {input}
        cd {params.dir}
        ($PASAHOME/scripts/Load_Current_Gene_Annotations.dbi -c {params.config} -g {params.genome} -P {params.anno2update}) > {log}
        cd ../../../../../
        touch {output}
        """

rule pasa_compare_anno1:
    input:
        loading="output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/load1.done",
        cfg="output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/annotCompare.done"
    output:
        "output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/compare1.done"
    log:
        config['workdir'] + "log/pasa/pasa_compare_anno1/{v}/{sample}.{run}.{tool}.pasa_{params}.pasa_compare_anno1.{v}.log"
    params:
        dir="output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/",
        config="annotCompare.config",
        genome=config['workdir'] + "/" + config['ref_genomes'],
        fasta_clean=config['workdir'] + "/output/isoform_generation/collapse_{run}/{v}/{sample}.collapsed.{run}.{v}.fasta.clean"
    threads: 15
    priority: 50
    conda:      
        "envs/pasa.yml"
    benchmark:
        "stats/pasa/pasa_compare_anno1/{v}/{sample}.{run}.{tool}.pasa_{params}.pasa_compare_anno1.{v}.txt"
    shell:
        """
        touch {input.loading}
        touch {input.cfg}
        cd {params.dir}
        ($PASAHOME/Launch_PASA_pipeline.pl -c {params.config} -A -g {params.genome} -t {params.fasta_clean} --CPU {threads}) > {log}
        cd ../../../../../
        touch {output}
        """

# rule pasa_load_anno2:
#     input:
#         check="output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/compare1.done",
#         config="output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/alignAssembly.config",
#         genome="input/refs/{sample}.scaffolds-v2.3.fasta"
#     output:
#         "output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/alignment2.done"
#     log:
#         "log/pasa/{v}/{run}/{sample}_pasa_{params}_loadAnno2.log"
#     params:
#         anno2update="input/annotations/{tool}/{sample}_mydb_pasa.gene_structures_post_PASA_updates.*.gff3"
#     priority: 40
#     shell:
#         "touch {input.check} && ~/software/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi -c {input.config} -g {input.genome} -P {params.anno2update} |& tee {log}"
# 
# rule pasa_compare_anno2:
#     input:
#         check="output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/alignment2.done",
#         config="output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/annotCompare.config",
#         genome="input/refs/{sample}.scaffolds-v2.3.fasta",
#         fasta_clean="output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/{sample}.collapsed.{run}.{v}.fasta.clean"
#     output:
#         "../../../../output/pasa/{run}/{v}/{sample}_{tool}_pasa_{params}/compare2.done"
#     log:
#         "log/pasa/{v}/{run}/{sample}_pasa_{params}_compAnno2.log"
#     params:
#         threads=30
#     priority: 30
#     shell:
#         "touch {input.check} && ~/software/PASApipeline/Launch_PASA_pipeline.pl -c {input.config} -A -g {input.genome} -t {input.fasta_clean} --CPU {params.threads} |& tee {log} && touch {output}"
