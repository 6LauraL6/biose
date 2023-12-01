from Bio.Seq import Seq

seq = Seq("GATCGATGGGCCCTATATAGGATCGAAAATCGC")

assert len(seq) == 33
assert seq.count("G") == 9

gc_fraction = (seq.count("G") + seq.count("C")) / len(seq)
assert round(gc_fraction, 2) == 0.48
