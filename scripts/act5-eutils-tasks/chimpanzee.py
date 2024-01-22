from Bio import Entrez
import json
from time import time
from tqdm import tqdm
from urllib.error import HTTPError



term = "chimpanzee[ORG] AND biomol mrna[PROP]"
db = "nucleotide"
file = "data/chimpanzee_mrna.fasta"
batch_size = 500

Entrez.email = "david@xtec.dev"

with Entrez.esearch(db,term, retmode = "json", usehistory = "y") as response:
  search = json.loads(response.read())["esearchresult"]


count = int(search["count"])
with open(file, "w") as writer:
    for start in tqdm(range(0, count, batch_size), desc=f"Fetch {count} (BS={batch_size}) ..."):
        end = min(count, start+batch_size)
       
        attempts = 0
        while attempts < 3:
            try:
                with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", retstart=start,        retmax=batch_size, webenv=search["webenv"], query_key=search["querykey"]) as response:
                    data = response.read()
                    writer.write(data)
                    break
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    attempts += 1
                    time.sleep(15)
                else:
                    raise