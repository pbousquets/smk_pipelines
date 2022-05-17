rule samtools_sort:
    input: 
        "{patient}/{sample}/{readgroup}/tmp.bam",
    output:
        temp("{patient}/{sample}/{readgroup}/tmp2.bam")
    resources:
        mem_mb=config.get("memory"),
    threads: 
        config.get("sort_threads", 4),
    singularity:
        "docker://labxa/bwamem2",
    shell:
        "samtools sort -@ {threads} -l 0 -m {resources.mem_mb}M -o {output} {input} 2>> .logs/sort.log"
