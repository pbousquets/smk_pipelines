rule BaseRecalibrator:
    input:
        bam =  "{sample}.dups.indelrealigned.bam",
        dbsnp = dbsnp,
        idx = idx
    output:
        "{sample}.recal_data.table",
    singularity:
        "docker://broadinstitute/gatk",
    shell:
        "gatk BaseRecalibrator -I {input.bam} -R {input.idx} --known-sites {dbsnp} -O {output}  2>> .logs/bqsr.log"

rule ApplyBQSR:
    input:
        bam =  "{sample}.dups.indelrealigned.bam",
        recal_table = "{sample}.recal_data.table",
        dbsnp = dbsnp,
        idx = idx
    output:
        "{sample}.dups.indelrealigned.bqsr.bam"
    singularity:
        "docker://broadinstitute/gatk",
    shell:
        "gatk --java-options '-Dsamjdk.compression_level=1' ApplyBQSR -R {idx} -I {input.bam} --bqsr-recal-file {input.recal_table} -O {output} 2>> .logs/bqsr.log"