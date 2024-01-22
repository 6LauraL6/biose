import json
import urllib.parse
import urllib3

db = "pubmed"
term = "science[journal] AND breast cancer AND 2008[pdat]"

#####

url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

term = urllib.parse.quote_plus(term)

response = urllib3.request(
    "GET",
    url,
    fields={"db": db, "term": term, "retmode": "json"}
)

data = response.json()
print(json.dumps(data["esearchresult"], indent=2))
