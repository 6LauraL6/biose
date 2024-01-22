from Bio import Entrez
import json

db = "pubmed"
term = "python and bioinformatics"
email = "david@xtec.dev"


def with_xml():

    record = Entrez.read(Entrez.esearch(db, term, email=email))

    for id in record["IdList"]:

        summary = Entrez.read(Entrez.esummary(db=db, id=id, email=email))
        summary = summary[0]
        print("Title = {}".format(summary["Title"]))
        print("DOI   = {}".format(summary["DOI"]))
        print("==============================================")

        quit()


def with_json():

    record = json.loads(Entrez.esearch(db, term, email=email, retmode="json").read())

    for id in record["esearchresult"]["idlist"]:

        summary = json.loads(Entrez.esummary(db=db, id=id, email=email, retmode="json").read())
        summary = summary["result"][id]
        print("Title = {}".format(summary["title"]))
        print("DOI   = {}".format(summary["elocationid"]))
        print("==============================================")

        quit()


with_json()
