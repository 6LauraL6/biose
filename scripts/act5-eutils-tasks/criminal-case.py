from Bio import Entrez, SeqIO
import os

gb_accession_prefix = 'AY' # First genbank is AY156734
gb_accession_ids = list(range(156734, 156764, 1))
save_folder = "./data/criminal-case"
Entrez.email = "mamoro10@xtec.cat"

print("Step 1. ePost Accession ID's.")
# Convert id numbers strings y agregar el prefijo
accession_ids_str = [f"{gb_accession_prefix}{id}" for id in gb_accession_ids]

search_handle = Entrez.epost(db="nucleotide", id=",".join(accession_ids_str))
search_results = Entrez.read(search_handle)
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]

print("Step 2. eFetch sequences info.")
handle = Entrez.efetch(db="nucleotide", webenv=webenv, query_key=query_key, rettype="gb", retmode="text")
genbank_records = list(SeqIO.parse(handle, "genbank"))

for i, record in enumerate(genbank_records):
    filename = os.path.join(save_folder, f"{gb_accession_prefix}{gb_accession_ids[i]}.gb")
    print(filename)
    #if not os.path.exists(filename):
    SeqIO.write(record, filename, "genbank")

print('Genbank files downloaded in ', save_folder)
