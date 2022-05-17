checkpoint checkpoint:
    input:
        bam=expand("{keys}/sorted.dups.bam",keys=paths)
    output:
        clusters=temp(touch("donedups"))
