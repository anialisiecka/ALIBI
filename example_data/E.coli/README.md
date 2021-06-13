You can use [this script](./runBio.py) to download genome records from NCBI using Biopython's [Bio.Entrez package](https://biopython.org/docs/1.75/api/Bio.Entrez.html).

It is assumed that you have a text file of NCBI GenBank accession numbers, one accession number per line (see [acc_numbers.txt](./acc_numbers.txt) as an example).

```
python3 runBio.py your_email_address text_file_of_accession_numbers
```
