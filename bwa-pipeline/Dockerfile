FROM broadinstitute/picard:2.27.1
RUN cd / && \
    apt update && \
    apt-get install curl bzip2 build-essential libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev -y && \
    curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2   --output - | tar -jxf -  && \
    curl -L https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2  --output - | tar -jxf -  && \
    cd /samtools-1.15.1  && \
    ./configure  && \
    make  && \
    make install   

RUN apt-get install python3 python3-pip -y && \
	pip3 install snakemake

RUN pip3 install pandas 

ENV PATH="${PATH}:/bwa-mem2-2.2.1_x64-linux" 
WORKDIR /usr/working
