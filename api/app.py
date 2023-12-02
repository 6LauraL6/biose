from flask import Flask
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
