FROM ubuntu

#Assemblers: spades, metaspades, megahit
#Gene Prediction:glimmer
#blast and BBMAP

#Needs viral dn nonviral dtabases
#location /home/cputonti/virmineAndrea/databases/
#all_viral_faa.fasta no_phage_bact_complete.fasta
RUN apt-get update && apt-get install -y python3 python3-dev gcc g++ unzip make git bzip2 zlib1g-dev ncurses-dev wget python3-pip ipython3 build-essential python3-pkg-resources python3-setuptools ncbi-blast+
#ADD all_viral_faa.fasta all_viral_faa.fasta
#ADD no_phage_bact_complete.fasta no_phage_bact_complete.fasta
ADD glimmer302b.tar.gz glimmer
ADD BBMap_37.36.tar.gz bbmap
ADD virMine.py virMine.py 
ADD SPAdes-3.10.1-Linux.tar.gz spades

#attempt to install biopython
RUN python3 -m pip install biopython


RUN git clone https://github.com/voutcn/megahit.git
#RUN tar xzf spades.tgz && tar xzf glimmer.tgz && tar xzf bbmap.tgz

#RUN mv spades/SPAdes-3.10.1-Linux/* spades/
#RUN rm -r spades/SPAdes-3.10.1-Linux

RUN pip3 install --upgrade pip


RUN cd megahit && make

ENV PATH /spades/bin:/megahit:/blast:$PATH
CMD ["python3", "virMine.py", "-h"]
RUN which python3
#sudo docker run -i -t thatzopoulos/virmine`
