rule RealignerTargetCreator:
    input:
        bam = "{patient}.merged.bam",
        dbsnp = dbsnp,
        idx = idx
    output:
        temp("{patient}.intervals")
    message:
        "Realigning indels!",
    threads: 
        config.get("samtools_threads", 4),
    singularity:
        "docker://broadinstitute/gatk3:3.8-1",
    shell:
        "java -Xmx24g -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator -known {input.dbsnp} -I {input.bam} -o {output} -R {input.idx} -nt {threads}"

rule GenomeAnalysisTK:
    input:
        bam = "{patient}.merged.bam",
        intervals = "{patient}.intervals",
        idx = idx
    output:
        "{patient}.realigned.bam"
    message:
        "Realigning indels!",
    threads: 
        config.get("samtools_threads", 4),
    singularity:
        "docker://broadinstitute/gatk3:3.8-1",
    shell:
        "java -Xmx24g -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner -I {input.bam}  -targetIntervals {input.intervals} -o {output} -R {input.idx} -compress 0"