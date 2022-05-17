def aggregate_input(wildcards):
    checkpoint_output = checkpoints.checkpoint.get(**wildcards).output[0]
    return checkpoint_output
    
rule samtools_merge:
    input:
        bams="{patient}",
        touch=aggregate_input
    output:
        temp("{patient}.merged.bam")
    threads: 
        config.get("other_threads", 4),
    singularity:
        "docker://labxa/bwamem2",
    shell:
        "samtools merge -@ {threads} -l 0 {output} {input.bams}/*/*/sorted.dups.bam && samtools index -@ {threads} {output} 2>> .logs/samtools_merge.log"
