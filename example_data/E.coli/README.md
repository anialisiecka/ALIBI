###### Step 1: Download genome records form NCBI
You can use [this script](./runBio.py) to download genome records from NCBI using Biopython's [Bio.Entrez package](https://biopython.org/docs/1.75/api/Bio.Entrez.html).

It is assumed that you have a text file of NCBI GenBank accession numbers, one accession number per line (see [acc_numbers.txt](./acc_numbers.txt) as an example). To run the script execute the following shell command:

```
python3 runBio.py <your_email_address> <text_file_of_accession_numbers>
```
The corresponding fasta file will be downloaded. The file will carry the same name as the input file, but with .fa extension.

###### Step 2: Convert FASTA file to GFA file
You can use the [vg toolkit](https://github.com/vgteam/vg) to construct sequence graph from your fasta file.

```
# make sequence graph
/path/to/vg/executable msga -f file_name.fa > file_name.vg

# convert the graph into GFA format
/path/to/vg/executable view file_name.vg > file_name.gfa
```
