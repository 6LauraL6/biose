import sys
from Bio.Seq import Seq
from Bio import SeqIO

file = "sequences.fa"
if len(sys.argv) != 2:
    print("Instruccions ús " + sys.argv[0] + " <fitxer fasta>")
    print(f"Com que no n'has proporcionat cap consultarem aquest: {file}")
else:
    file = sys.argv[1]

# Llegim fitxer
sequences = SeqIO.parse(open(file), "fasta")
print(f'Llegint el fitxer amb aminoàcids {file}')

# Obtenim le seqüències i les transformem a proteïnes
for record in sequences:
   print("Seq Prot: ", record.id, record.seq)