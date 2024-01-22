import urllib.parse
import urllib3

db = "pubmed"
term = "science[journal] AND breast cancer AND 2008[pdat]"

#####

# https://www.ncbi.nlm.nih.gov/books/NBK25500/
url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

term = urllib.parse.quote_plus(term)
response = urllib3.request(
    "GET",
    url,
    fields={"db": db, "term": term})

data = response.data.decode("utf-8")
print(data)
