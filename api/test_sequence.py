from Bio.Seq import Seq

def test_dna_rna():
   coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

   template_dna = coding_dna.reverse_complement()
   messenger_rna = template_dna.reverse_complement().transcribe()

   assert messenger_rna == Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")

# def test2():
#     seq = Seq("GATCGA")
#     assert seq.reverse_complement() == Seq("TCGATT")