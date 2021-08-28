# ALIBI - Algorithm for Linearization by Incremental graph BuIlding

## Prerequisites
Running:
* [Python 3](https://www.python.org/)
* [Biopython](https://biopython.org/)

## Installation
1. Clone this repository:
```
git clone https://github.com/anialisiecka/ALIBI.git
```
2. Change directory:
```
cd ALIBI
```
3. Run the setup.py file from that directory:
```
python setup.py install
```

## Input
ALIBI requires as input a GFA file (Graphical Fragment Assembly) or a MAF file (Multiple Alignment File) with a sequence graph.

## Output
ALIBI linearizes the genome sequence graph present in the input file. The final linearized graph is written to a file of the same format as the input file.

## Usage
To use ALIBI you need to be in the ALIBI folder.
```
bash alibi.sh -i <input_file>
```
The final linearized graph is written to a GFA or MAF file with the same name as the input file, but with "\_sorted" suffix before the file extension.

## Founding
This software is developed with support of [OPUS 11 scientific project of National Science Centre: Incorporating genomic variation information into DNA sequencing data analysis.](https://www.mimuw.edu.pl/~dojer/rmg/)

## License
This project is licensed under the MIT License - see the [LICENSE](./LICENSE) file for details.
