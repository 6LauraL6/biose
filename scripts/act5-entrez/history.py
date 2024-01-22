from Bio import Entrez
from tqdm import tqdm

term = "Opuntia[orgn] and rpl16"
file = "data/orchid_rpl16.fasta"
batch_size = 3

Entrez.email = "david@xtec.dev"

with Entrez.esearch("nuccore", term, usehistory="y") as response:
    record = Entrez.read(response)

count = int(record["Count"])
with open(file, "w") as writer:
    for start in tqdm(range(0, count, batch_size), desc=f"Fetch {count} (BS={batch_size}) ..."):
        end = min(count, start + batch_size)

        with Entrez.efetch(db="nuccore", rettype="fasta", retmode="text", retstart=start, retmax=batch_size, webenv=record["WebEnv"], query_key=record["QueryKey"]) as response:
            writer.write(response.read())
