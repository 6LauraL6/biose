from Bio                import SeqIO            # SeqIO is a module (several classes)
from Bio.SeqRecord      import SeqRecord        # SeqRecord is a class
#from Bio.SeqFeature     import SeqFeature, FeatureLocation # Inner classes
from pprint             import pformat, pp      # Així veiem millor els diccionaris.

# Llegim el fitxer, que és unigenbank.
GENBANK_NAME: str = "NC_005129.2_elephas_mitoc.genbank"
with open(GENBANK_NAME, "r") as file:
    demo_genbank: SeqRecord = SeqIO.read(GENBANK_NAME, "genbank")
    # També funciona amb:
    # demo_genbank: list[SeqRecord] = list(SeqIO.read(GENBANK_NAME, "genbank"))
    # Resum del genbank.
    # print(demo_genbank)

    print('SeqRecord from genbank file:')
    # Com hem fet amb el fasta, podem llegir dades del SeqRecord del genbank.
    # Annotacions, del bloc 1. 
    # Ens retorna un diccionari. Els camps compostos (com autors) són una llista.
    pp(demo_genbank.annotations)

    # Si només ens interessen les referències i els seus autors.
    authors = []
    for reference in demo_genbank.annotations.get('references', []):
        authors.extend(reference.authors.split(", "))
    # Si volem eliminar duplicats
    authors = set(authors)
    print("Llista d'autors:", authors)

    # Solució més completa:
    references = demo_genbank.annotations.get('references', [])
    for i, reference in enumerate(references, 1):
        authors = reference.authors.split(", ")
        title = reference.title
        journal = reference.journal
        pubmed = reference.pubmed_id

        print(f"Referència {i}:")
        print(f"  Autors: {authors}")
        print(f"  Títol: {title}")
        print(f"  Journal: {journal}")
        print(f"  PubMed ID: {pubmed}")

    # Si volem obtenir totes les features en general.
    # Per no allargar la sortida, posem les 3 primeres.
    features = demo_genbank.features[0:3]
    for feature in features:
        print("Tipo de característica:", feature.type)
        print("Ubicació:", feature.location)

        # Accedeix als qualifiers (etiquetes) de la característica
        qualifiers = feature.qualifiers
        for key, value in qualifiers.items():
            print(f"{key}: {value}")

        print()
        
    # Si només ens interessen els comentaris.
    comments = demo_genbank.annotations.get('comment', None)
    if comments:
        print("Llista de comentaris:", comments)

    # Per llegir seqüència de nucleòtids, del bloc 3, és igual que el fasta.
    # Fixeu-vos que surt tota junta, a diferència del genbank.
    print('Primers 100 caràcters Seqüència ')
    print(demo_genbank.seq[0:100])
