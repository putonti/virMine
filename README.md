# virMine

Novel viral species genomes – particularly those in high abundance – have been able to be identified directly from complex community viral metagenomes. Discovery of such viral genomes often relies heavily on manual curation and prior studies have employed a variety of different criteria when sifting through sequencing data. In an effort to provide a comprehensive means for the discovery of complete viral genomes from complex sequence data sets, we developed the tool virMine. Raw sequencing reads are processed and assembled and annotated, and individual contigs are scored based upon their likelihood of being viral in origin. Several filters have been implemented, allowing researchers to refine their search for their specific study system. virMine can be used to identify viruses in any niche and thus further our understanding of this vast reservoir of genetic diversity.

## Getting Started

Clone Project

```python
git clone https://github.com/thatzopoulos/virmine.git
```

Move paired-end fastq files, as well as the viral and nonviral databases, to the inputFiles folder prior to building the docker image

From within the project folder run:
```python
sudo docker build -t virmine .
```
```python
sudo docker run -v /pathToLocalFolder/runName:/virmineDockerOutputFolder -i -t virmine
```

### Prerequisites

Docker is the only prerequisite for this program to run, all other dependencies are handled by the Dockerfile.

## virMine Command Options

*	a : choose your assembler (spades, metaspades, megahit)
*	p : list your paired-end read files
*	v : list your viral database
*	nv : list your nonviral database
*	o : make an output file

### Example Run 
```python
python2.7 virMine.py -a spades -p inputFiles/R1.fastq inputFiles/R2.fastq -v inputFiles/viral_aa.fasta -nv inputFiles/nonviral_aa.fasta -o outputFolder
```

## Authors

* Andrea Garretto
* Thomas Hatzopoulos
* Catherine Putonti

## License

This project is licensed under the terms of the MIT License
