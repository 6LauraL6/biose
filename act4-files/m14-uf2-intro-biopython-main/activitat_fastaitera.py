from Bio import SeqIO                           # SeqIO is a module
from Bio.SeqIO.FastaIO  import FastaIterator    # FastaIterator is a class

ORCHID_FASTA = 'ls_orchid.fasta'

# When not enough RAM, use iterators to read one SeqRecord 
# at a time (memory-efficient)
record_iterator: FastaIterator = SeqIO.parse(ORCHID_FASTA, 'fasta')

# Iterators can be traversed only once.
# All SeqRecords have at least the .seq attribute (Seq object)
print('Lengths of all orchids:')
for record in record_iterator:
    print( len(record.seq), end=', ' )
