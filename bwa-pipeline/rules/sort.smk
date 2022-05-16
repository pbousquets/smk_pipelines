rule samtools_sort:
    input: 
        "{patient}/{sample}/{readgroup}/tmp.bam",
    output:
        temp("{patient}/{sample}/{readgroup}/tmp2.bam")
    log:
        ".logs/{patient}/{sample}/{readgroup}/sort.log",
    message:
        "Sorting!",
    resources:
        mem_mb=100000,
    threads: 
        config.get("samtools_threads", 4),
    singularity:
        "docker://labxa/bwamem2",
    shell:
        "samtools sort -@ {threads} -l 0 -m {resources.mem_mb}M -o {output} {input} 2>{log}"
