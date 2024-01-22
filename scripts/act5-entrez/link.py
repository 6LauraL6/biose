from Bio import Entrez
from tabulate import tabulate

pmid = "14907713"

Entrez.email = "david@xtc.dev"

with Entrez.esummary(db="pubmed", id=pmid) as response:
  links = Entrez.read(response)[0]
  #print(record.keys())
  print(f'{links["LastAuthor"]} ({links["PubDate"]}). {links["Title"]} {links["ELocationID"]}')
  print("")

with Entrez.elink(dbfrom="pubmed", id=pmid) as response:
    links = Entrez.read(response)[0]

links = [[link["DbTo"], link["LinkName"], len(link["Link"])] for link in links["LinkSetDb"]]

print(tabulate(links, headers=["DbTo", "LinkName", "len(Link)"]))


