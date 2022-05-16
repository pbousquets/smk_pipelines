rule mergeRG:
    input:
        touch = "donesplit",
    output:
        "{sample}.dups.indelrealigned.bam",
    threads: 
        config.get("samtools_threads", 7),
    singularity:
        "docker://labxa/bwamem2",
    shell:
        """
        nfiles=$(ls {wildcards.sample}.*.dups.indelrealigned.bam | wc -l)
        if [ $nfiles -gt 1 ]
        then
            samtools merge -@ {threads} -l 0 {output} {wildcards.sample}.*.dups.indelrealigned.bam
        else
            mv {wildcards.sample}.*.dups.indelrealigned.bam {output}
        fi && samtools index -@ {threads} {output}
        """
