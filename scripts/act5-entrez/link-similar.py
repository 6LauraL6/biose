from Bio import Entrez

pmid = "14907713"

Entrez.email = "david@xtc.dev"

def print_summary(summary):
    print(
        f'{summary["LastAuthor"]} ({summary["PubDate"]}). {summary["Title"]} {summary["ELocationID"]}')

with Entrez.esummary(db="pubmed", id=pmid) as response:
  summary = Entrez.read(response)[0]
  print_summary(summary)


with Entrez.elink(dbfrom="pubmed", id=pmid) as response:
    record = Entrez.read(response)[0]

similar = [ls for ls in record["LinkSetDb"]
           if ls["LinkName"] == "pubmed_pubmed"][0]
ids = [link["Id"] for link in similar["Link"]]

with Entrez.esummary(db="pubmed", id=",".join(ids)) as response:
    record = Entrez.read(response)

print("\nSimilar articles\n=============================================================")

for summary in record[0:3]:
    print_summary(summary)


