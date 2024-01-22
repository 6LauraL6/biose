from Bio import Entrez, SeqIO
import os.path

def buscar_accession_id(gen, organismes):
    accession_ids = {}

    for organisme in organismes:
        # Realitzar una cerca per l'organisme i el gen específic
        query = f'{organisme}[Organism] AND {gen}[Gene]'
        handle = Entrez.esearch(db='nucleotide', term=query, retmax=1)
        result = Entrez.read(handle)

        # Afegir l'Accession ID al diccionari
        if 'IdList' in result and result['IdList']:
            accession_ids[organisme] = result['IdList'][0]

    return accession_ids

def descarregar_genbank_si_cal(accession_ids):
    for organisme, accession_id in accession_ids.items():
        # Verificar si el fitxer ja existeix localment
        filename = f'{organisme.replace(" ", "_")}_{accession_id}.gb'
        if os.path.isfile(filename):
            print(f"El fitxer GenBank per {organisme} ja existeix: {filename}")
        else:
            # Descarregar el fitxer GenBank utilitzant l'Accession ID
            handle = Entrez.efetch(db='nucleotide', id=accession_id, rettype='gb', retmode='text')
            genbank_data = handle.read()

            # Guardar el fitxer GenBank amb un nom específic
            with open(filename, 'w') as file:
                file.write(genbank_data)
            print(f"Fitxer GenBank per {organisme}: {filename}")


# Configurar les API keys per accedir als serveis de NCBI
Entrez.email = "mamoro10@xtec.cat"

# Organismes d'interès i gen comú (p53)
organismes_interes = [
    "Rattus norvegicus",
    "macaca fascicularis",
    "Canis lupus familiaris"
]

gen_comu = "p53"

# Buscar Accession IDs
accession_ids = buscar_accession_id(gen_comu, organismes_interes)

# Descarregar fitxers GenBank
descarregar_genbank_si_cal(accession_ids)
