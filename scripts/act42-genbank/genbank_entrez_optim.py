from Bio import Entrez
from Bio import SeqIO
import os

# Substitueix amb el teu nom de fitxer, base de dades NCBI i Accession.
GENBANK_NAME = "SarsCov2.gb"  
DATABASE = "nuccore"
ACCESSION = "NC_045512"

# GENBANK_NAME = "HomoSapiens.gb"  
# DATABASE = "protein"
# ACCESSION = "4507341"

empremtes_guardades = os.path.join(os.getcwd(), GENBANK_NAME)
Entrez.email = "mamoro10@xtec.cat"

if not os.path.exists(empremtes_guardades):
    # Descarrega el genbank si no existeix
    #with Entrez.efetch(db="protein", id='4507341', rettype="gb", retmode="text") as handle:
    with Entrez.efetch(db=DATABASE, id=ACCESSION, rettype="gb", retmode="text") as handle:
        genbank_content = handle.read()
    
    # Guarda el contingut al fitxer genbank
    with open(empremtes_guardades, 'w') as genbank_file:
        genbank_file.write(genbank_content)
    
    print(f"Genbank descarregat i guardat com a {GENBANK_NAME}")
else:
    print(f"Ja existeix un fitxer Genbank amb el nom {GENBANK_NAME}. Procedim a llegir-lo.")
    # Llegeix multifasta o multigenbank, retorna un Iterador.
    demo_genbank = list(SeqIO.parse(GENBANK_NAME, "genbank"))
    print(demo_genbank)
