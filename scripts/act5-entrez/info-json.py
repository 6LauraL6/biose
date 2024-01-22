from Bio import Entrez
import json
from tabulate import tabulate

Entrez.email = "david@xtec.dev"

db = "nuccore"

info = json.loads(Entrez.einfo(db=db, retmode="json").read())
fieldlist = info["einforesult"]["dbinfo"][0]["fieldlist"]

fields = [[field["name"], field["fullname"], field["termcount"],
           field["description"]] for field in fieldlist]

print(tabulate(fields, headers=[
      "name", "fullname", "termcount",  "description"]))

'''
{
  "name": "ALL",
  "fullname": "All Fields",
  "description": "All terms from all searchable fields",
  "termcount": "13471894110",
  "isdate": "N",
  "isnumerical": "N",
  "singletoken": "N",
  "hierarchy": "N",
  "ishidden": "N"
}
'''
