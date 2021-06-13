import sys
import os
from Bio import Entrez, SeqIO

arguments = sys.argv

Entrez.email = arguments[1]
file_name = arguments[2]

input_file = open(file_name,"r")
output_file = open(os.path.splitext(file_name)[0]+".fa","w")
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
