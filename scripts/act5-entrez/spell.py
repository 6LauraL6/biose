from Bio import Entrez

term = "adenin"

Entrez.email = "david@xtec.dev"

with Entrez.espell(term = term) as response:
  record = Entrez.read(response)

print(record["CorrectedQuery"])


