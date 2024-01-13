from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import os

# Directori on es troben els fitxers fasta i genbank
FASTA_DIR = os.path.join(os.path.dirname(__file__), 'data')
GENBANK_DIR = os.path.join(os.path.dirname(__file__), 'data')

# Session 2 tools.

def processar_seq_adn(seq_adn):
    # Verificar que les bases siguin correctes
    seq_valida = all(base in "ACGTacgt" for base in seq_adn)

    if seq_valida:
        # Normalitzar majúscules i minúscules
        seq_normalitzada = Seq(seq_adn.upper())

        # Obtindre l'ARNm
        arnm = seq_normalitzada.transcribe()

        # Obtindre l'ADN plantilla
        plantilla_adn = seq_normalitzada.reverse_complement()

        # Calcular el percentatge GC
        percentatge_gc = GC(seq_normalitzada)

        return {
            "success": True,
            "adn_plantilla": str(plantilla_adn),
            "arnm": str(arnm),
            "percentatge_gc": percentatge_gc
        }
    else:
        return {
            "success": False,
            "error": "La seqüència d'ADN conté bases no vàlides."
        }
    
# Session 3 tools.

def translate_mrna(mrna_sequence, translation_table=1):
    mrna_seq = Seq(mrna_sequence)
    protein_seq = mrna_seq.translate(table=translation_table)
    return protein_seq

def analyze_cds(protein_seq):
    start_codon = "M"  # El codón de inicio es representado por "M" en la secuencia de aminoácidos
    stop_codons = ["*", "Ter"]  # Los codones de parada son representados por "*" o "Ter" en la secuencia de aminoácidos

    cds_segments = []
    incomplete_proteins = []

    for start_index in range(0, len(protein_seq)):
        amino_acid = protein_seq[start_index]
        if amino_acid == start_codon:
            cds = protein_seq[start_index:]
            for stop_codon in stop_codons:
                stop_index = cds.find(stop_codon)
                if stop_index != -1:
                    cds_segments.append(str(cds[:stop_index + 1]))
                    break
            else:
                incomplete_proteins.append(str(cds))

    return cds_segments, incomplete_proteins

# Session 4 tools.

def get_fasta_info(code):
    fasta_file = f'{FASTA_DIR}/{code}.fasta'
    
    try:
        records = list(SeqIO.parse(fasta_file, 'fasta'))
        print(records[0].seq[0:30])  # Debug: Imprimeix la seqüència del primer registre
        print(f'Reading {fasta_file}')
        info = {
            'descriptions': [record.description for record in records],
            'sequences': [str(record.seq) for record in records]
        }
        return info
    except FileNotFoundError:
        print(f'File not found: {fasta_file}')  # Debug: Imprimeix la ruta completa del fitxer
        return {'error': f'No es troba cap fitxer fasta amb el codi {code}. Suggereixo provar amb NC_045512 (el del SarsCov2).'}


def get_genbank_files_list():
    genbank_files = [file for file in os.listdir(GENBANK_DIR) if file.endswith('.genbank')]
    return genbank_files


def get_genbank_info(accession):
    genbank_file = f'{GENBANK_DIR}/{accession}'
    print(genbank_file)
    try:
        records = SeqIO.parse(genbank_file, 'genbank')
        info = []
        
        for record in records:
            record_info = {
                'title': record.description,
                'accession': record.id,
                'ncbi_link': f'https://www.ncbi.nlm.nih.gov/nuccore/{record.id}',
                'latest_reference': record.annotations.get('references', [])[0].title,
                'num_features': len(record.features),
                'cds_info': [{'type': feature.type, 'location': str(feature.location)} for feature in record.features if feature.type == 'CDS'],
                'sequence': str(record.seq),
                'gc': round(GC(record.seq))
            }
            info.append(record_info)

        return info
    except FileNotFoundError:
        return {'error': f'No es troba cap fitxer genbank amb l\'accession {accession}. Suggerim provar amb el NC_045512; el del SarsCov2.'}