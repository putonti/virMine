# virMine

TODO

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

* TODO
* Thomas Hatzopoulos
* Catherine Putonti

## License

This project is licensed under the TODO

## Acknowledgments

* TODO
