from Bio                import SeqIO            # SeqIO is a module
from Bio.Seq            import Seq              # Seq is a class
from Bio.SeqRecord      import SeqRecord        # SeqRecord is a class
from Bio.SeqFeature     import SeqFeature       # Class to manage genbank features.

from pathlib import Path        # Manage file and folder paths.
from pprint  import pp         # Show dicts clearer way.

# Constants.
DATA_DIR = Path('./data')

# Fasta and Genbank source, last checked in 08/01/2024:
# https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2
SARS_COV_2_FASTA    = DATA_DIR/'sars-cov-2.fasta'
SARS_COV_2_GENBANK  = DATA_DIR/'sars-cov-2.genbank'

# #############################################################################
# - Biopython Exercises
# #############################################################################
# q0 - Mètodes per llegir el genbank.
# q1 - Quina longitud té la secuencia ? Mostra’n els 30 primers caràcters.
# q2 - Mostrar quantes annotations del context (bloc 1) hi ha, 
#      l’accession i la referència més recent (la primera del fitxer genbank).
# q3 - Mostra el type i la location de totes les features del genbank.
# q4 - Analitza la última seqüència i compara la seqüència del genbank amb la del fasta.

# -----------------------------------------------------------------------------
def q0_read_SarsCov2_genbank() -> SeqRecord:
    '''Read and return default Genbank as a SeqRecord.'''
    sars_cov_2_record: SeqRecord = SeqIO.read(SARS_COV_2_GENBANK, 'genbank')
    return sars_cov_2_record

def q0_read_other_genbank(genbank_name: str) -> SeqRecord:
    '''Read and return a genbank file with name defined in genbank_name as a SeqRecord.'''
    '''genbank_name: str -> Name of the genbank to read.'''
    sars_cov_2_record: SeqRecord = SeqIO.read(genbank_name, 'genbank')
    return sars_cov_2_record
# -----------------------------------------------------------------------------
    
def q1(genbank_rec: SeqRecord) -> str:
    '''Return the sequence. Show seq. length and first 30 nucleotides.'''

    print("TODO")
    return "CAGACAGAT"

# -----------------------------------------------------------------------------
def q2(genbank_rec: SeqRecord):
    '''Show which annotations are in the genbank file, the accession, and the last reference'''

    print("TODO")

# -----------------------------------------------------------------------------
def q3(genbank_rec: SeqRecord):
    '''Show how many features are in the genbank file.'''

    print("TODO")


# -----------------------------------------------------------------------------
def q4(genbank_rec: SeqRecord):
    '''Get the last CDS, translate it manually and compare it to the genbank translation.'''

    print("TODO")
    #print(f"Same translation? {same_translation}")


# -----------------------------------------------------------------------------
def main():

    sars_cov_2_genbank: SeqRecord = q0_read_SarsCov2_genbank()
    # No estem obligats a agafar el paràmetre de retorn de la q1 (str: sequencia)
    q1(sars_cov_2_genbank)
    q2(sars_cov_2_genbank)
    q3(sars_cov_2_genbank)
    q4(sars_cov_2_genbank)


# Main
# -----------------------------------------------------------------------------
this_module = __name__
main_module = '__main__'

if (this_module == main_module): main()
# -----------------------------------------------------------------------------