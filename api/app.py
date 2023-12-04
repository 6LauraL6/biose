from flask import Flask, request, jsonify
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import sequence
import time

app = Flask(__name__, static_folder="out", static_url_path="/")


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




if __name__ == "__main__":
    app.run()


from Bio import SeqIO
from Bio.Seq import Seq

seq = Seq("AGTA")

assert Seq("TCAT") == seq.complement()


def load_seq():
    # curl https://www.ebi.ac.uk/ena/browser/api/fasta/AP006852.1?download=true -o ap006852.fasta
    records = list(SeqIO.parse("ap006852.fasta", "fasta"))
    dna = records[0]

    print(dna.name)
    print(dna.description)
    print(dna.seq[:100])
