rule RealignerTargetCreator:
    input:
        bam = "{patient}.merged.bam",
        dbsnp = dbsnp,
        idx = idx
    output:
        temp("{patient}.intervals")
    threads: 
        config.get("other_threads", 4),
    singularity:
        "docker://broadinstitute/gatk3:3.8-1",
    shell:
        "java -Xmx24g -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator -known {input.dbsnp} -I {input.bam} -o {output} -R {input.idx} -nt {threads} 2>> .logs/indel_realigner.log"

rule GenomeAnalysisTK:
    input:
        bam = "{patient}.merged.bam",
        intervals = "{patient}.intervals",
        idx = idx
    output:
        "{patient}.realigned.bam"
    singularity:
        "docker://broadinstitute/gatk3:3.8-1",
    shell:
        "java -Xmx24g -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner -I {input.bam}  -targetIntervals {input.intervals} -o {output} -R {input.idx} -compress 0 2>> .logs/indel_realigner.log"