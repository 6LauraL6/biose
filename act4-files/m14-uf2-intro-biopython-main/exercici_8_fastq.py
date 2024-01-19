import sys
from Bio import SeqIO

# Calcular probabilitats error cada caràcter del fastq
def calcular_probabilitats_error(puntuacions_qualitat):
    probabilitats = [10**(-ord(caracter) / 10.0) for caracter in puntuacions_qualitat]
    return probabilitats

file = "sequences.fq"
# Llegeix un fastq per defecte si no li hem posat cap per paràmetre.
if len(sys.argv) != 2:
    print("Instruccions d'ús: " + sys.argv[0] + " <fitxer fasta d'ADN o ARN>")
    print(f"Com que no n'has proporcionat cap, consultarem aquest: {file}")
else:
    file = sys.argv[1]

# llegeix el nom del fitxer FASTQ
data = SeqIO.parse(open(file), "fastq")
print(f'Reading file {file}')
print("File "+file+" readed succesfully.")

# Calculem puntuacions qualitat.
puntuacions_qualitat = ['!', 'A', 'a', '%']
probabilitats_error = calcular_probabilitats_error(puntuacions_qualitat)

for i in range(len(puntuacions_qualitat)):
    print(f"Puntuació: {puntuacions_qualitat[i]}, Probabilitat d'error: {probabilitats_error[i]}")
