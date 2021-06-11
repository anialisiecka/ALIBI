You can use the follownig code to download a list of gene IDs from NCBI.

```python
from Bio import Entrez, SeqIO
Entrez.email="test@test.org"

input_file = open("example.txt","r")
output_file = open("output.fa","w")
acc = ""

for line in input_file:
        line = line.strip()
        acc += line
        acc += ","
acc = acc[:-1]

input_file.close()

handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
records = SeqIO.parse(handle, "fasta")

SeqIO.write(records,output_file,"fasta")

handle.close()
output_file.close()
```
