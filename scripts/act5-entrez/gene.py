from Bio import Entrez

db = "gene"
term = "cobalamin synthase homo sapiens"
email = "david@xtec.dev"

with Entrez.esearch(db, term, email=email) as response:
  record = Entrez.read(response)

for id in record["IdList"]:

    with Entrez.esummary(db=db, id=id, email=email) as response:
      summary = Entrez.read(response)
    
    summary = summary["DocumentSummarySet"]["DocumentSummary"][0]
    #print(summary.keys())

    print("Id          = {}".format(id))
    print("Name        = {}".format(summary["Name"]))
    print("Description = {}".format(summary["Description"]))
    print("Summary     = {}".format(summary["Summary"]))
    print("==============================================")

    quit()
