include: "utils.smk"
rule bwa_mem2:
    input:
        unpack(get_fastq),
        idx = idx
    output:
        pipe("{patient}/{sample}/{readgroup}/tmp.bam"),
    log:
        ".logs/{patient}/{sample}/{readgroup}/bwa.log"
    params:
        rg = r"-R '@RG\tID:{sample}.{readgroup}\tSM:{patient}\tPL:ILLUMINA\tPU:ILLUMINA\tLB:WGS'", #Keep sample in RG ID to allow the pipeline merge back later
    message:
        "Aligning!"
    threads: 
        config.get("bwa_threads", 4),
    singularity:
        "docker://labxa/bwamem2",
    shell:
        "bwa-mem2 mem -t {threads} {params.rg} {input.idx} {input.fastq1} {input.fastq2} > {output} 2>{log}"