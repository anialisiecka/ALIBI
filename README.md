# ALIBI - Algorithm for Linearization by Incremental graph BuIlding

## Input
ALIBI requires as input a GFA file with a sequence graph.

## Output
ALIBI linearizes the genome sequence graph present in the input file. The final linearized graph is written to a GFA file.

## Usage
Let ```/path/to/alibi/``` denote the path to ALIBI root directory. To run ALIBI execute the following shell command: 
```
bash /path/to/alibi/bin/alibi.sh -i <input_file>
```
The final linearized graph is written to a GFA file with the same name as the input file, but with "\_sorted" suffix before the file extension.

## Founding
This software is developed with support of [OPUS 11 scientific project of National Science Centre: Incorporating genomic variation information into DNA sequencing data analysis.](https://www.mimuw.edu.pl/~dojer/rmg/)

## License
This project is licensed under the MIT License - see the [LICENSE](./LICENSE) file for details.
