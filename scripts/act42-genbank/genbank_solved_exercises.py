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

    print(f'q1 -> Sequence', len(genbank_rec.seq))
    print(f'Number of gene features:', genbank_rec.seq[0:30])
    return genbank_rec.seq

# -----------------------------------------------------------------------------
def q2(genbank_rec: SeqRecord):
    '''Show which annotations are in the genbank file, the accession, and the last reference'''

    print(f'q2 -> Number of context annotations:', len(genbank_rec.annotations))
    # Debug: pp(f'Annotations:', genbank_rec.annotation)
    print(f'Accession ID:', genbank_rec.annotations['accessions'][0])
    # Last reference
    print(f'Last reference data: \n', genbank_rec.annotations['references'][0])

# -----------------------------------------------------------------------------
def q3(genbank_rec: SeqRecord):
    '''Show how many features are in the genbank file.'''

    print(f"q3 -> Number of features: {len(genbank_rec.features)}")

    # Show type and location of each feature
    feature_info_list: list[SeqFeature] = [ (feature.type, str(feature.location))
                                            for feature
                                            in genbank_rec.features ]

    pp(feature_info_list, width=200)


# -----------------------------------------------------------------------------
def q4(genbank_rec: SeqRecord):
    '''Get the last CDS, translate it manually and compare it to the genbank translation.'''

    # Get list of CDS features
    cds_list: list[SeqFeature] = [  feature
                                    for feature
                                    in genbank_rec.features
                                    if feature.type == 'CDS']

    # Get last cds
    cds_feature: SeqFeature = cds_list[-1]

    # Explore the CDS SeqFeature
    pp(cds_feature)

    # Get start and end
    cds_start: int = int(cds_feature.location.start)
    cds_end:   int = int(cds_feature.location.end)

    # Get my translation with BioPython
    cds_seq = genbank_rec.seq[cds_start:cds_end]
    my_translation: Seq = cds_seq.translate()

    # Get genbank translation
    genbank_translation: Seq = Seq(cds_feature.qualifiers['translation'][0])

    # Checks
    print(my_translation)
    print(genbank_translation)

    # Genbank translation with stop
    genbank_translation_with_stop: Seq = genbank_translation + '*'

    # Same translation?
    same_translation = (my_translation == genbank_translation_with_stop)
    print(f"Same translation? {same_translation}")


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