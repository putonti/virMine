FROM ubuntu

#Assemblers: spades, metaspades, megahit
#Gene Prediction:glimmer
#blast and BBMAP

#Location: /home/cputonti/virmineAndrea/databases
#Viral DB: all_viral_faa.fasta
#NonViral DB: no_phage_bact_complete.fasta

RUN apt-get update && apt-get install -y curl vim csh python3 python3-dev gcc g++ unzip make git bzip2 zlib1g-dev ncurses-dev wget python3-pip ipython3 build-essential python3-pkg-resources python3-setuptools ncbi-blast+
#ADD glimmer302b.tar.gz glimmer
ADD BBMap_37.36.tar.gz bbmap
ADD virMine.py virMine.py 
ADD SPAdes-3.10.1-Linux.tar.gz spades
ADD inputFiles/ /inputFiles/
#attempt to install biopython
RUN pip3 install --upgrade pip
RUN python3 -m pip install biopython


RUN git clone https://github.com/voutcn/megahit.git
RUN git clone https://github.com/najoshi/sickle
#RUN tar xzf glimmer.tgz && tar xzf bbmap.tgz

#RUN mv spades/SPAdes-3.10.1-Linux/* spades/
#RUN rm -r spades/SPAdes-3.10.1-Linux/

RUN cd sickle && make
RUN cd megahit && make

#RUN cd glimmer/glimmer3.02/src && make

#GET Glimmer
ENV GLIMMER_VERSION 302b
ENV GLIMMER_DIR /opt/glimmer
ENV GLIMMER_SUBDIR glimmer3.02
RUN mkdir -p $GLIMMER_DIR
RUN curl -SL https://ccb.jhu.edu/software/glimmer/glimmer$GLIMMER_VERSION.tar.gz | tar -xzC $GLIMMER_DIR
RUN cd $GLIMMER_DIR/$GLIMMER_SUBDIR/src && make
#Add glimmer to PATH
ENV PATH $GLIMMER_DIR/$GLIMMER_SUBDIR/bin:$GLIMMER_DIR/$GLIMMER_SUBDIR/scripts:$PATH
#update paths in the g3-iterated.csh script
#set awkpath = /glimmer/glimmer3.02/scripts
#set glimmerpath = /glimmer/glimmer3.02/bin

#Get ELPH
ENV ELPH_VERSION 1.0.1
ENV ELPH_DIR /opt/ELPH
RUN mkdir -p $ELPH_DIR
RUN curl -SL ftp://ftp.cbcb.umd.edu/pub/software/elph/ELPH-$ELPH_VERSION.tar.gz | tar -xzC $ELPH_DIR
RUN cd $ELPH_DIR/ELPH/sources && make
RUN mkdir -p $ELPH_DIR/bin
RUN mv $ELPH_DIR/ELPH/sources/elph $ELPH_DIR/bin/elph
ENV PATH $ELPH_DIR/bin:$PATH
#Update hard-coded paths in g3-iterated.csh 
RUN sed -i "s|/fs/szgenefinding/Glimmer3|${GLIMMER_DIR}/${GLIMMER_SUBDIR}|g" $GLIMMER_DIR/$GLIMMER_SUBDIR/scripts/g3-iterated.csh
RUN sed -i "s|/nfshomes/adelcher/bin/elph|${ELPH_DIR}/bin/elph|g" $GLIMMER_DIR/$GLIMMER_SUBDIR/scripts/g3-iterated.csh
#update awk paths in awk scripts
RUN sed -i "s|/bin/awk|/usr/bin/awk|g" $GLIMMER_DIR/$GLIMMER_SUBDIR/scripts/*.awk
#add bbmap to path
ENV PATH /bbmap/bbmap:/spades/SPAdes-3.10.1-Linux/bin:/megahit:/blast:/sickle:$PATH
RUN echo $PATH
#CMD ["python3", "virMine.py", "-a", "spades", "-p", "inputFiles/R1.fastq", "inputFiles/R2.fastq", "-v", "inputFiles/viral_aa.fasta", "-nv", "inputFiles/nonviral_aa.fasta", "-o", "testOutput"]
RUN which python3
#sudo docker run -i -t thatzopoulos/virmine`
