from Bio import Entrez
import json

email = "david@xtec.dev"


def with_xml():

    # Return a list of all Entrez database names:
    with Entrez.einfo(email=email) as response:
        record = Entrez.read(response)
        assert ("pubmed" in record["DbList"])

    # Return statistics for Entrez PubMed:
    with Entrez.einfo(db="pubmed", email=email) as response:
        record = Entrez.read(response)["DbInfo"]
        assert (record["DbName"] == "pubmed")
        print(record["Count"])

def with_json():
    
    with Entrez.einfo(retmode="json", email=email) as response:
        record = json.load(response)
        assert ("pubmed" in record["einforesult"]["dblist"])

    with Entrez.einfo(db = "pubmed", retmode="json", email=email) as response:
        record = json.load(response)
        #print(json.dumps(record, indent= 2))

        info = record["einforesult"]["dbinfo"][0]
        assert(info["dbname"] == "pubmed")
        print(info["count"])

def list_fields():
    
  with Entrez.einfo(db="pubmed", email=email) as response:
        record = Entrez.read(response)["DbInfo"]
        
        print("Name\tFullName\tDescription")
        print("-----------------------------------------------------------------------")
        for field in record["FieldList"]:
          print("%(Name)s\t'%(FullName)s'\t %(Description)s" % field)

list_fields()
