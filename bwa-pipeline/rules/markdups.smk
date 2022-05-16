rule picard_markdups:
    input:
        "{patient}/{sample}/{readgroup}/tmp2.bam",
    output:
        bam="{patient}/{sample}/{readgroup}/sorted.dups.bam",
        metrics="{patient}.{sample}.{readgroup}.marked_dup_metrics.txt",
    log:
        ".logs/{patient}/{sample}/{readgroup}/picard_dups.log",
    message:
        "Marking dups!",
    threads: 1,
    singularity:
        "docker://broadinstitute/picard",
    shell:
        "java -jar /usr/picard/picard.jar MarkDuplicates -I {input} -O {output.bam} -M {output.metrics}  --COMPRESSION_LEVEL 0 2>{log}"