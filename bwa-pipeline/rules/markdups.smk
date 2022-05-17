rule picard_markdups:
    input:
        "{patient}/{sample}/{readgroup}/tmp2.bam",
    output:
        bam="{patient}/{sample}/{readgroup}/sorted.dups.bam",
        metrics="{patient}.{sample}.{readgroup}.marked_dup_metrics.txt",
    singularity:
        "docker://broadinstitute/picard",
    shell:
        "java -jar /usr/picard/picard.jar MarkDuplicates -I {input} -O {output.bam} -M {output.metrics}  --COMPRESSION_LEVEL 0 2>> .logs/picard_dups.log"