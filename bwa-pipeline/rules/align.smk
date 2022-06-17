include: "utils.smk"
rule bwa_mem2:
    input:
        unpack(get_fastq),
        idx = idx
    output:
        pipe("{patient}/{sample}/{readgroup}/tmp.bam"),
    params:
        PL = config['PL'],
        PU = config['PU'],
        CN = config['CN'],
        rg = r"-R '@RG\tID:{sample}.{readgroup}\tSM:{patient}\tPL:{params.PL}\tPU:{params.PU}\tLB:{sample}\tCN:{params.CN}'", #Keep sample in RG ID to allow the pipeline merge back later
    threads: 
        config.get("bwa_threads", 18),
    singularity:
        "docker://labxa/bwamem2",
    shell:
        "bwa-mem2 mem -t {threads} {params.rg} {input.idx} {input.fastq1} {input.fastq2} > {output} 2>> .logs/bwa.log"
