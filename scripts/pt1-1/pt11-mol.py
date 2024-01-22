import os
from urllib.request import urlretrieve

folder_url = "./data/"

# URL Wikipedia.
# format: "https://es.wikipedia.org/wiki/${nom}"
url_wikipedia_en: str  = "https://en.wikipedia.org/wiki/"
url_wiki_completed_en: str  = ""

# Paràmetre URL, ha de ser un nom de macromolecula.
# Nosaltres volem obtenir les següents:
adn_mol_names: list[str] = ["Adenine","Cytosine","Thymine","Guanine"]
# adn_mol_urls: list[str]  = []

print("Obtaining PUBchem macromolecules")

# ... your code ...

print("PUBchem macromolecules obtained successfully!")