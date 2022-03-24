# virMine 2.0

Novel viral species genomes – particularly those in high abundance – have been able to be identified directly from complex community viral metagenomes. Discovery of such viral genomes often relies heavily on manual curation and prior studies have employed a variety of different criteria when sifting through sequencing data. In an effort to provide a comprehensive means for the discovery of complete viral genomes from complex sequence data sets, we developed the tool virMine. Input can be either: (1) raw sequencing reads from short-read sequencing technologies or (2) long-read sequences or assembled contigs. Short-reads are processed and assembled. These assembled short-reads or supplied long-reads/contigs are then annotated, and individual contigs are scored based upon their likelihood of being viral in origin. Several filters have been implemented, allowing researchers to refine their search for their specific study system. virMine can be used to identify viruses in any niche and thus further our understanding of this vast reservoir of genetic diversity.

## Getting Started

Clone Project

```python
git clone https://github.com/putonti/virmine.git
```

Move paired-end fastq files, as well as the viral and nonviral databases, to the inputFiles folder prior to building the docker image

From within the project folder run:
```python
sudo docker build --tag virmine:latest virmine
```
```python
sudo docker run -v ~/pathToLocalFolder/virmine:/virmineDockerOutputFolder -i -t virmine
```

### Prerequisites

Docker is the only prerequisite for this program to run, all other dependencies are handled by the Dockerfile.
If any section of the program causes an error or is unable to run, check that you have enough memory in your Docker resources.

## virMine Command Options

*	-a : choose your assembler (spades, metaspades, megahit, all3)
*	-p : list your short read files (one single-end file or two paired-end files) 
*	-v : list your viral database
*	-nv : list your nonviral database
*	-o : make an output file

If you wish to supply long-reads or assembled contigs:
*	-A : list your long-read/contig assembly file

### Example Run with Paired-End Reads:
```python
python3 virMine.py -a spades -p inputFiles/R1.fastq inputFiles/R2.fastq -v inputFiles/viral_aa.fasta -nv inputFiles/nonviral_aa.fasta -o virmineDockerOutputFolder/output
```
### Example Run with Long-reads or Assembled Contigs:
```python
python3 virMine.py -A inputFiles/assembled_contigs.fasta -v inputFiles/viral_aa.fasta -nv inputFiles/nonviral_aa.fasta -o virmineDockerOutputFolder/output
```

## Test Data
Artificial paired-end reads of _E. coli_ K12 (NC_000913.3) and Escherichia phage Lambda (J02459.1).

### Example Command to Run Test Data:
```python
python3 virMine.py -a all3 -p testFiles/test_data_R1.fastq testFiles/test_data_R2.fastq -v inputFiles/viral_aa.fasta -nv inputFiles/nonviral_aa.fasta -o virmineDockerOutputFolder/testOutput
```

### Example Output
Running the test data using the option of all 3 assembly methods should result in 3 contigs (final_contigs.fasta), with 1 of those contigs being classified as viral (viral_contigs.fasta). Other output files include all predicted ORFs, the viral ORFs only, and all BLAST output files.  

## Further Details:
Further details about this tool can be found in Garretto et al. PeerJ 2019 (https://peerj.com/articles/6695/). A manuscript describing the upgrades implemented as part of v2.0 is currently under review.

## Authors

* Andrea Garretto
* Thomas Hatzopoulos
* Genevieve Johnson
* Catherine Putonti

## License

This project is licensed under the terms of the MIT License
