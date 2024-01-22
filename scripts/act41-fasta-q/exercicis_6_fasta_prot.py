import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

file = "sequences.fa"
if len(sys.argv) != 2:
    print("Instruccions d'ús: " + sys.argv[0] + " <fitxer fasta d'ADN o ARN>")
    print(f"Com que no n'has proporcionat cap, consultarem aquest: {file}")
else:
    file = sys.argv[1]

# Llegeix fitxer
try:
    sequences = SeqIO.parse(open(file), "fasta")
    print(f'Llegint el fitxer {file}')
except Exception as e:
    print(f"No s'ha pogut llegir el fitxer: {e}")
    sys.exit(1)

# Inicialitza una llista per emmagatzemar SeqRecords de proteïnes
seq_records_proteina = []

# Recorre cada seq_record i tradueix la seqüència a proteïna
for seq_record in sequences:
    try:
        seq_proteina = seq_record.seq.translate()
        print(seq_proteina)
        nova_capcalera = f"{seq_record.description} translated by Biopython"
        print(nova_capcalera)
        seq_record_proteina = SeqRecord(seq_proteina, id=seq_record.id, description=nova_capcalera)
        seq_records_proteina.append(seq_record_proteina)
    except Exception as e:
        print(f"Excepció en processar {seq_record.id}: {e}")
        sys.exit(1)

# Guarda els resultats
output_file = file.replace(".fasta", "_newprot.fasta")
SeqIO.write(seq_records_proteina, output_file, "fasta")
print(f"Les seqüències traduïdes s'han guardat a {output_file}")