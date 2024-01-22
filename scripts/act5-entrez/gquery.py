from Bio import Entrez
from tabulate import tabulate

Entrez.email = "david@xtec.dev"

with Entrez.egquery(term="SARS-CoV-2") as response:
  record = Entrez.read(response)

result = [[row["DbName"], row["Count"]] for row in record["eGQueryResult"] ]
print(tabulate(result, headers=["DbName","Count"]))
