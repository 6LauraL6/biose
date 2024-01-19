from Bio import Entrez

Entrez.email = "mamoro10@xtec.cat"
handle = Entrez.efetch(db="protein", id='4507341', rettype="gb", retmode="text")
print(handle.read())