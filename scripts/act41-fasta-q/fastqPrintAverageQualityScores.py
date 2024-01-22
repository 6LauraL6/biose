import sys
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as pyplot

# aquesta secció s'encarrega de llegir les dades de l'usuari
if len(sys.argv) != 2 or "-h" in sys.argv or "--help" in sys.argv:
    print("Usage: fastqPrintAverageQualityScores.py <fastq_file>")
    sys.exit()

# llegeix el nom del fitxer FASTQ
fastq = sys.argv[1]
# analitza el fitxer FASTQ
data = SeqIO.parse(fastq, "fastq")

print("File "+fastq+" readed succesfully.")

print("File "+fastq+" content:")
for record in data:
   print(record.description)
   print(record.seq)

print("Calculating quality.")
# inicialitza una llista de probabilitats
sum_p = [0] * 200
N = [0] * 200

for record in data:
    for i, Q in enumerate(record.letter_annotations["phred_quality"]):
        # converteix la puntuació PHRED a una probabilitat
        p_err = 10**(-float(Q) / 10.0)
        # afegeix aquesta probabilitat específica a l'array per a aquesta seqüència
        sum_p[i] += p_err
        N[i] += 1

# ara per a la representació gràfica. Inicialitza les llistes x i y per a la representació gràfica:
x = []
y = []

# afegeix les mitjanes de les probabilitats als valors y:
for i in range(len(N)):
    if N[i] > 0:
        pAvg = sum_p[i] / N[i]
        x.append(i)
        y.append(pAvg)

# representa gràficament els valors x i y:
print("Generating quality plot in file quality.png")
pyplot.plot(x, y)
pyplot.xlabel('posició (nt)')
pyplot.ylabel('Probabilitat mitjana d\'error')
pyplot.savefig('quality.png')

print("Plot generated successfully.")
