checkpoint checkpoint:
    input:
        bam=expand("{keys}/sorted.dups.bam",keys=paths)
    output:
        clusters=touch("donedups")
