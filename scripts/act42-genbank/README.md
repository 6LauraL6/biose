## Lectura de Biopython.

La teoria que proporcionem respecte Biopython es troba en aquest document (esborrany):

https://docs.google.com/document/d/1j5ayfY9b3tz7tFeErp7EJkWqnTohJ0VHvu6kNo7rTkg/edit

És la continuació d'aquest repositori d'introducció a Biopython:

https://github.com/miquelamorosaldev/m14-uf2-intro-biopython/

### Altres Referències:

- https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec33
- https://www.ncbi.nlm.nih.gov/gene

## Exemples de com llegir fitxers de Genbank.

Passem a resumir-ho amb codi.

### Es pot fer manualment.

[read_genbank_manual.py](https://github.com/miquelamorosaldev/m14-uf2-intro-genbank/blob/main/read_genbank_manual.py)

Genbank de prova:
[NC_005129.2_elephas_mitoc.genbank](https://github.com/miquelamorosaldev/m14-uf2-introbiopython-genbank/blob/main/NC_005129.2_elephas_mitoc.genbank)

### Amb la llibreria Entrez.

Solució senzilla:
[read_genbank_entrez.py](https://github.com/miquelamorosaldev/m14-uf2-intro-genbank/blob/main/read_genbank_entrez.py)

Solució més robusta:
[genbank_entrez_optim.py](https://github.com/miquelamorosaldev/m14-uf2-intro-genbank/blob/main/genbank_entrez_optim.py)

Si volem descarregar-nos diversos recursos a partir d'una llista d'AccessionId aquí tenim exemples de com fer-ho:
[sum-eutils.py](https://github.com/miquelamorosaldev/m14-uf2-intro-genbank/blob/main/sum-eutils.py)

Es combinen les Entrez eUtils d'aquesta manera: 
* eSearch  -> Obtenir llista Accession
* ePost    -> Publiquem llista Accession i obtenin sessionId
* eFetch   -> Usem el sessionId per descarregar-nos els fitxers.

## Exercicis.

La resta de fitxers són exercicis bàsics per practicar tot el que es pot llegir del Genbank amb BioPython; enunciats i solucions.

## Resum contingut fitxers Genbank.

### Bloc 1 - Context (o Annotations)
* 29903 posició de la seqüencia en la que acaba la seqüència.
* DBLINK projecte d'on ha sortit aquesta seqüència. 
* ORGANISM Diferents classificacions del coronavirus. 
* Reference Authors tots els que l'han seqüenciat 
* Title Títol del article on han publicat. 
* JOURNAL Mitjà on ho han publicat 
* PUBMED link del article (si el volem llegir podem intentar buscar a SCI HUB) 
* REFERENCE de otras personas que han estat treballant sobre el tema. 
* FEATURES les diferents anotacions que tenen el fitxer. [anotacions]

El fitxer // es el final d'un fitxer genbank

### Bloc 2 - Features (anotacions dels gens)

En tot cas, descrivim breument els camps que té dins del Coronavirus:

* source.	Dades generals sobre la cadena: longitud, organisme, tipus…
* gene		Informació del gen
* CDS 		Cadena codificant /translation Cadena proteïna. Podem comprovar que s’ha codificat bona part del virus.
* Sovint hi ha molts CDS. A vegades s’han fet joins de diversos CDS.
Dins de cada CDS tenim info útil com el gen, el codó d’inici (si és l’1 es tradueix amb la M), l’id de la proteïna funcional.

Podem trobar altres features:
* 5'UTR 	En aquest cas veiem que des de la 1..265 no ho ha codificat gene.
* mat_peptide 	Regions concretes de la proteïna que és important remarcar, ja que fan funcions dins la seqüència (enllaç amb altres proteïnes per exemple).


### Bloc 3 - ORIGIN. Seqüència nucleòtids completa

// -> fi del fitxer


