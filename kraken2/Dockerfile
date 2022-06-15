from staphb/kraken2

ENV TZ=Europe/Madrid
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN \
	apt-get update && \
	apt-get install -y libbz2-dev liblzma-dev && \
	cd / && wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2 && \
	tar xvjf samtools-1.15.1.tar.bz2 && \
	cd samtools-1.15.1 && \
	./configure --without-curses && make && make install && \
	cd /


RUN \ 
	mkdir /db && \
	cd /db && \
	wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20220607.tar.gz && \
	tar -xf k2_viral_20220607.tar.gz && \
	rm k2_viral_20220607.tar.gz && \
	cd /

RUN apt-get install python3 python3-pip -y && \
        pip3 install snakemake

RUN	cd / && apt-get install -y git && \
	git clone https://github.com/jenniferlu717/Bracken  && \
	cd Bracken && bash install_bracken.sh 

RUN pip3 install pandas click

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 0

ENV PATH="/Bracken:${PATH}"

COPY . /kraken_pipeline/

ENTRYPOINT ["/kraken_pipeline/wrapper.py"]
