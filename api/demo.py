from Bio.Seq import Seq

seq = Seq("GATCGA")
result = str(seq.reverse_complement())
print("Seqüència Original:", seq)
print("Invers Complementari de la seqüència:", result)
