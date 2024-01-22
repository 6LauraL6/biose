from Bio import Entrez
from time import time
from tqdm import tqdm
from urllib.error import HTTPError

term = "Opuntia[ORGN]"
file = "data/orchid_papers.txt"
batch_size = 10

Entrez.email = "david@xtec.dev"

with Entrez.esearch(db="pubmed", term=term, reldate=365, datetype="pdat", usehistory="y") as response:
    search = Entrez.read(response)

count = int(search["Count"])
with open(file, "w") as writer:
    for start in tqdm(range(0, count, batch_size), desc=f"Fetch {count} (BS={batch_size}) ..."):
        end = min(count, start+batch_size)
       
        attempts = 0
        while attempts < 3:
            try:
                with Entrez.efetch(db="pubmed", rettype="medline", retmode="text", retstart=start,        retmax=batch_size, webenv=search["WebEnv"], query_key=search["QueryKey"]) as response:
                    data = response.read()
                    writer.write(data)
                    break
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    attempts += 1
                    time.sleep(15)
                else:
                    raise
