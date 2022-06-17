rule bwa_mem2:
    input:
        unpack(get_fastq),
        idx = idx
    output:
        pipe("{readgroup}_tmp.bam"),
    params:
        PL = config['PL'],
        SM = config['SM'],
        LB= config['LB'],
        CN = config['CN'],
    threads: 
        config.get("bwa_threads", 18),
    singularity:
        "docker://labxa/bwamem2",
    log:
        ".logs/{readgroup}/bwa.log"
    shell:
        """
        PU=$(zcat ${input.fastq1} | head -1 | sed 's/[:].*//' | sed 's/@//') && \\
        rg = r"-R '@RG\tID:{readgroup}\tSM:{params.SM}\tPL:{params.PL}\tPU:{PU}\tLB:{params.LB}\tCN:{params.CN}'", #Keep sample in RG ID to allow the pipeline merge back later
        bwa-mem2 mem -t {threads} {params.rg} {input.idx} {input.fastq1} {input.fastq2} > {output} 2>> {log}"""
