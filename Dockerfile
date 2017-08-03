FROM ubuntu

#Assemblers: spades, metaspades, megahit
#Gene Prediction:glimmer
#blast and BBMAP

#Location: /home/cputonti/virmineAndrea/databases
#Viral DB: all_viral_faa.fasta
#NonViral DB: no_phage_bact_complete.fasta

RUN apt-get update && apt-get install -y python3 python3-dev gcc g++ unzip make git bzip2 zlib1g-dev ncurses-dev wget python3-pip ipython3 build-essential python3-pkg-resources python3-setuptools ncbi-blast+
<<<<<<< HEAD
=======
ADD all_viral_faa.fasta all_viral_faa.fasta
#ADD no_phage_bact_complete.fasta no_phage_bact_complete.fasta
>>>>>>> runs, needs input files
ADD glimmer302b.tar.gz glimmer
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
<<<<<<< HEAD
#RUN rm -r spades/SPAdes-3.10.1-Linux/
=======
#RUN rm -r spades/SPAdes-3.10.1-Linux

>>>>>>> runs, needs input files

RUN cd sickle && make
RUN cd megahit && make
<<<<<<< HEAD
RUN cd glimmer/glimmer3.02/src && make
#add bbmap to path
ENV PATH /bbmap/bbmap:/glimmer/glimmer3.02/bin:/spades/SPAdes-3.10.1-Linux/bin:/megahit:/blast:/sickle:$PATH
RUN echo $PATH
#CMD ["python3", "virMine.py", "-a", "spades", "-p", "inputFiles/R1.fastq", "inputFiles/R2.fastq", "-v", "inputFiles/viral_aa.fasta", "-nv", "inputFiles/nonviral_aa.fasta", "-o", "testOutput"]
=======
#add bbmap to path
ENV PATH /spades/bin:/megahit:/blast:$PATH
CMD ["python3", "virMine.py", "-h"]
>>>>>>> runs, needs input files
RUN which python3
#sudo docker run -i -t thatzopoulos/virmine`
