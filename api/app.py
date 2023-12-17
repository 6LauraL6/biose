from flask import Flask, request, jsonify, Response, json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC, seq3
import sequence
import time

app = Flask(__name__, static_folder="out", static_url_path="/")

## Web Services Sessió 2.

@app.route("/")
def index():
    return app.send_static_file("index.html")


@app.route("/api/time")
def get_current_time():
    return {"time": time.time()}

@app.route("/api/sequence/<seq>")
def get_sequence_info(seq):
    info = sequence.info(seq)

    return {"info": info}

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

@app.route('/processar_adn', methods=['POST'])
def processar_adn():
    data = request.get_json()

    if 'seq_adn' in data:
        seq_adn = data['seq_adn']
        resultat = processar_seq_adn(seq_adn)
        return jsonify(resultat)
    else:
        return jsonify({"success": False, "error": "La seqüència d'ADN no s'ha proporcionat."})

## Web Services Sessió 3.

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

@app.route('/translate', methods=['POST'])
def translate():
    data = request.get_json()

    mrna_sequence = data.get('mrna_sequence', '')
    # Verificar que les bases siguin correctes
    arn_seq_valid = all(base in "ACGUacgu" for base in mrna_sequence)

    if arn_seq_valid:
        # Normalitzar majúscules i minúscules
        seq_normalitzada = Seq(mrna_sequence.upper())
        # Taula traducció estandard, es pot enviar una altra.
        translation_table = data.get('translation_table', 1)

        protein_seq = translate_mrna(mrna_sequence, translation_table)
        cds_segments, incomplete_proteins = analyze_cds(protein_seq)

        amino_acids = {acid: seq3(acid) for acid in str(protein_seq)}

        result = {
            'success': True,
            'protein_sequence': str(protein_seq),
            'cds_segments': cds_segments,
            'incomplete_proteins': incomplete_proteins,
            'amino_acids': amino_acids,
        }
    else:
        result = {
            "success": False,
            "error": "La seqüència d'ADN conté bases no vàlides."
        }
    json_string = json.dumps(result,ensure_ascii = False)
    response = Response(json_string,content_type="application/json; charset=utf-8" )
    return response

## Web Services necessaris per arrencar la aplicació.

if __name__ == "__main__":
    app.run()

# seq = Seq("AGTA")
# assert Seq("TCAT") == seq.complement()

def load_seq():
    # curl https://www.ebi.ac.uk/ena/browser/api/fasta/AP006852.1?download=true -o ap006852.fasta
    records = list(SeqIO.parse("ap006852.fasta", "fasta"))
    dna = records[0]

    print(dna.name)
    print(dna.description)
    print(dna.seq[:100])
