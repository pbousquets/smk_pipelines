include: "utils.smk"
rule bwa_mem2:
    input:
        unpack(get_fastq),
        idx = idx
    output:
        pipe("{patient}/{sample}/{readgroup}/tmp.bam"),
    params:
        PL = config['PL'],
        SM = config['SM'],
        LB= config['LB'],
        CN = config['CN'],
        ID = r"{readgroup}", #Keep sample in RG ID to allow the pipeline merge back later
    threads: 
        config.get("bwa_threads", 18),
    singularity:
        "docker://labxa/bwamem2",
    shell:
        """
        pu=$(zcat {input.fastq1} | head -1 | sed 's/[:].*//' | sed 's/@//') && \ 
        rg="@RG\\tID:{params.ID}\\tSM:{params.SM}\\tPL:{params.PL}\\tPU:$pu\\tLB:{params.LB}\\tCN:{params.CN}" && \
        bwa-mem2 mem -t {threads} -R $rg {input.idx} {input.fastq1} {input.fastq2} > {output} 2>> {log}
        """
