from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import sys

# Exercici 3 i 4. Llegeix fitxer fasta.
file = "sequences.fa"
if len(sys.argv) != 2:
    print("Instruccions d'ús: " + sys.argv[0] + " <fitxer fasta d'ADN o ARN>")
    print(f"Com que no n'has proporcionat cap, consultarem aquest: {file}")
else:
    file = sys.argv[1]

sequences = SeqIO.parse(open(file), "fasta")
print(f'Reading file {file}')
list_seq_code = []

# Exercici 5. Validació de bases.
# Si alguna no és correcte es sortirà del programa.
for record in sequences:
   #print(record.id, record.seq)
   list_seq_code.append(record.seq)
   seq_aux = Seq(record.seq)
   # Exercici 1.
   seq_rev = seq_aux.reverse_complement()
   print(f'Secuence = {seq_aux[:10]}')
   print(f'Reversed = {seq_rev[:10]}')
   # Exercici 2.
   gc = gc_fraction(seq_aux)*100
   print(f'GC = {gc} %')
   # Exercici 5.
   is_unambiguous = all(base in "ACGTacgt" for base in seq_aux)
   print(f'Is unambiguous? {is_unambiguous}')
   
   