rule bam_split:
    input:
        expand("{patient}.realigned.bam", patient=pats)
    output:
        touch("donesplit")
    threads: 
        config.get("other_threads", 4),
    singularity:
        "docker://labxa/bwamem2",
    shell:
        "samtools split -@ {threads} {input} -f \%\!.dups.indelrealigned.bam 2>> .logs/bam_split.log"

