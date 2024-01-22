from Bio import Entrez
from tqdm import tqdm

ids = list(range(19000000, 19002000))
file = "data/pubmed.txt"
batch_size = 500

Entrez.email = "david@xtec.dev"

id = ",".join(list(map(str, ids)))
with Entrez.epost("pubmed", id=id) as response:
    record = Entrez.read(response)

count = len(ids)
with open(file, "w") as writer:
    for start in tqdm(range(0, count, batch_size), desc=f"Fetch {count} (BS={batch_size}) ..."):
        end = min(count, start + batch_size)

        with Entrez.efetch(db="pubmed", rettype="medline", retmode="text", retstart=start, retmax=batch_size, webenv=record["WebEnv"], query_key=record["QueryKey"]) as response:
            writer.write(response.read())
