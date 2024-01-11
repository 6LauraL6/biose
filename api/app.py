from flask import Flask, request, jsonify, Response, json
import sequence
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
import time
import biotools
#from flask_cors import CORS

app = Flask(__name__, static_folder="out", static_url_path="/")
#cors = CORS(app, resources={r"/translate": {"origins": "http://localhost:3000"}})

@app.route("/")
def index():
    return app.send_static_file("index.html")

@app.route("/api/time")
def get_current_time():
    return {"time": time.time()}

## Web Services Sessió 2.

@app.route("/api/sequence/<seq>")
def get_sequence_info(seq):
    info = sequence.info(seq)
    return {"info": info}

@app.route('/api/processar_adn', methods=['POST'])
def processar_adn():
    data = request.get_json()

    if 'seq_adn' in data:
        seq_adn = data['seq_adn']
        resultat = biotools.processar_seq_adn(seq_adn)
        return jsonify(resultat)
    else:
        return jsonify({"success": False, "error": "La seqüència d'ADN no s'ha proporcionat."})

@app.route('/api/translate', methods=['POST'])
def translate():
    data = request.get_json()

    mrna_sequence = data.get('mrna_sequence', '')
    # Verificar que les bases siguin correctes
    arn_seq_valid = all(base in "ACGUacgu" for base in mrna_sequence)

    if arn_seq_valid:
        # Normalitzar majúscules i minúscules
        seq_normalitzada = Seq(mrna_sequence.upper())
        # Taula traducció estandard, es pot enviar una altra.
        # https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables
        translation_table = data.get('translation_table', 1)

        protein_seq = biotools.translate_mrna(mrna_sequence, translation_table)
        cds_segments, incomplete_proteins = biotools.analyze_cds(protein_seq)

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
    # response.headers.add('Access-Control-Allow-Origin', 'http://localhost:3000')  # Ajusta l'origen
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type')
    return response

## Web Services Sessió 4.

@app.route('/fasta/<code>', methods=['GET'])
def get_fasta(code):
    print(code)
    info = biotools.get_fasta_info(code)
    print(info)
    return jsonify(info)

@app.route('/genbank/<accession>', methods=['GET'])
def get_genbank(accession):
    info = biotools.get_genbank_info(accession)
    return jsonify(info)

@app.route('/api/genbank-list', methods=['GET'])
def get_genbank_list():
    genbank_files = biotools.get_genbank_files_list()
    genbank_info_list = []

    for accession in genbank_files:
        genbank_info = biotools.get_genbank_info(accession.split('.')[0])
        genbank_info_list.append({
            'title': genbank_info['title'],
            'accession': genbank_info['accession'],
            'ncbi_link': genbank_info['ncbi_link'],
            'latest_reference': genbank_info['latest_reference'],
            'num_features': genbank_info['num_features'],
            'seq_first30': genbank_info['seq_first30'],
            'cds_info': genbank_info['cds_info'],
            'repeated_sequence': genbank_info['repeated_sequence']
        })

    return jsonify(genbank_info_list)


## Web Services necessaris per arrencar la aplicació.

if __name__ == "__main__":
    app.run()

# seq = Seq("AGTA")
# assert Seq("TCAT") == seq.complement()

# def load_seq():
#     # curl https://www.ebi.ac.uk/ena/browser/api/fasta/AP006852.1?download=true -o ap006852.fasta
#     records = list(SeqIO.parse("ap006852.fasta", "fasta"))
#     dna = records[0]

#     print(dna.name)
#     print(dna.description)
#     print(dna.seq[:100])
