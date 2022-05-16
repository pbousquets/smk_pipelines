rule bam_split:
    input:
        expand("{patient}.realigned.bam", patient=pats)
    output:
        touch("donesplit")
    message:
        "Spliting!",
    threads: 
        config.get("samtools_threads", 7),
    singularity:
        "docker://labxa/bwamem2",
    shell:
        "samtools split -@ {threads} {input} -f \%\!.dups.indelrealigned.bam"

