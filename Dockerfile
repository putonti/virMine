FROM ubuntu

#Assemblers: spades, metaspades, megahit
#Gene Prediction:glimmer
#blast and BBMAP

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y curl vim csh python3.9 python3.9-dev gcc g++ unzip make git bzip2 zlib1g-dev ncurses-dev wget python3-pip build-essential python-pkg-resources python-setuptools ncbi-blast+ cmake
RUN pip install ipython 
ADD BBMap_38.94.tar.gz bbmap
ADD virMine.py virMine.py 
ADD SPAdes-3.15.3-Linux.tar.gz spades
ADD inputFiles/ /inputFiles/
ADD testFiles/ /testFiles/
ADD g3-iterated-viral.csh g3-iterated-viral.csh 
RUN pip install --upgrade pip
RUN python3 -m pip install biopython


RUN wget https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
RUN git clone https://github.com/najoshi/sickle


RUN cd sickle && make
RUN tar zvxf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz

#GET Glimmer
ENV GLIMMER_VERSION 302b
ENV GLIMMER_DIR /opt/glimmer
ENV GLIMMER_SUBDIR glimmer3.02
RUN mkdir -p $GLIMMER_DIR
RUN curl -SL --cipher DEFAULT@SECLEVEL=1 https://ccb.jhu.edu/software/glimmer/glimmer302b.tar.gz | tar -xzC $GLIMMER_DIR
RUN cd $GLIMMER_DIR/$GLIMMER_SUBDIR/src && make
#Add glimmer to PATH
ENV PATH $GLIMMER_DIR/$GLIMMER_SUBDIR/bin:$GLIMMER_DIR/$GLIMMER_SUBDIR/scripts:$PATH

#Get ELPH
ENV ELPH_VERSION 1.0.1
ENV ELPH_DIR /opt/ELPH
RUN mkdir -p $ELPH_DIR
RUN curl -SL --cipher DEFAULT@SECLEVEL=1 ftp://ftp.cbcb.umd.edu/pub/software/elph/ELPH-$ELPH_VERSION.tar.gz | tar -xzC $ELPH_DIR
RUN cd $ELPH_DIR/ELPH/sources && make
RUN mkdir -p $ELPH_DIR/bin
RUN mv $ELPH_DIR/ELPH/sources/elph $ELPH_DIR/bin/elph
ENV PATH $ELPH_DIR/bin:$PATH
RUN chmod u+x $ELPH_DIR/bin
RUN chmod u+x g3-iterated-viral.csh

#Update hard-coded paths in g3-iterated.csh 
RUN sed -i "s|/fs/szgenefinding/Glimmer3|${GLIMMER_DIR}/${GLIMMER_SUBDIR}|g" $GLIMMER_DIR/$GLIMMER_SUBDIR/scripts/g3-iterated.csh
RUN sed -i "s|/nfshomes/adelcher/bin/elph|${ELPH_DIR}/bin/elph|g" $GLIMMER_DIR/$GLIMMER_SUBDIR/scripts/g3-iterated.csh
RUN sed -i "s|/bin/awk|/usr/bin/awk|g" $GLIMMER_DIR/$GLIMMER_SUBDIR/scripts/*.awk

RUN cp $GLIMMER_DIR/$GLIMMER_SUBDIR/scripts/g3-iterated.csh .

ENV PATH /bbmap/bbmap:/spades/SPAdes-3.15.3-Linux/bin:/megahit:/blast:/sickle:$PATH
RUN echo $PATH

#Example Run Command
#CMD ["python3.9", "virMine.py", "-a", "spades", "-p", "inputFiles/R1.fastq", "inputFiles/R2.fastq", "-v", "inputFiles/viral_aa.fasta", "-nv", "inputFiles/nonviral_aa.fasta", "-o", "testOutput"]
