## Introducció a Biopython.

La teoria que proporcionem respecte Biopython es troba en aquest document (esborrany):

https://docs.google.com/document/d/1v0BcVDdad613peyC_aMyqM0ehR33idG-AhkUf1WbU8Y/edit

### Referències:

- https://open.oregonstate.education/appliedbioinformatics/chapter/chapter-1/
- https://biopython.org/


## Com llegir fitxers .fasta amb la shell de Linux:

Molt fàcil!

```bash
grep -E '^\S+$' sequences.fa | awk '{printf "%s ", $0; getline; print $0}'
```

## Exercicis bloc 1. Lectura fitxers fasta de nucleòtids.

[exercicis_1_5_fasta.py](https://github.com/miquelamorosaldev/m14-uf2-intro-biopython/blob/main/exercicis_1_5_fasta.py)

Escriu un programa que llegeixi un fitxer FASTA indicat per codi:

1.- Imprimeixi un fitxer FASTA nou amb el complement invers = reverse complementary de cada seqüència.

2.- El contingut GC de cada seqüència.

3.- Descarrega’t el fitxer .fasta de la orquidea, prova el codi i analitza els resultats.

https://github.com/biopython/biopython/blob/master/Doc/examples/ls_orchid.fasta

4.- Millora el programa (si encara no ho has fet) i fes que si l’usuari introdueix un fitxer fasta pel terminal 
que imprimeixi el fasta del terminal, i si no ho ha fet, que avisi a l’usuari que la idea és que posi un fitxer fasta, 
però en comptes de sortir del programa que llegeixi un fasta passat per codi (el secuences.fa)

5.- Afegeix un mecanisme per verificar quines seqüències del fitxer FASTA només tinguin seqüències d'ADN vàlides i sense ambigüitats: 
(validant els caràcters) i ho mostri per pantalla. 

<em>Observació: Tot i que ja vam veure que ens podem trobar seqüències amb caràcters com la N 
(que indica que pot ser qualsevol àcid nucleic = qualsevol lletra) i és un caràcter vàlid en els fitxers fasta; no ens interessa tenir cadenes ambigües.</em>

## Activitat. Lectura fitxers fasta d'aminoàcids.

[activitat_llegir_prot.py](https://github.com/miquelamorosaldev/m14-uf2-intro-biopython/blob/main/activitat_llegir_prot.py)

També és molt important recordar com llegir un fitxer .fasta de nucleotids i convertir-lo a proteïna; l‘objectiu de l’activitat és fer crear el codi que converteixi adn a proteïnes. 
Al final, aconsegueix que guardi el fitxer resultant en un altra fasta, que tingui el mateix nom excepte que els últims caràcters siguin “_newprot”.fasta

## Exercicis bloc 2. Lectura de nucleòtids i escriptura d'aminoàcids en firxers fasta.

[exercicis_6_fasta_prot.py](https://github.com/miquelamorosaldev/m14-uf2-intro-biopython/blob/main/exercicis_6_fasta_prot.py)

Seguint amb el pas 5 (solució de l'exercici), fes un nou programa que tradueixi un seqüència d’ADN o un ARN a una de proteïnes. 
Al final, ha de guardar la proteïna en un altre fitxer que tingui el mateix nom juntament amb el text “_newprot”. 
La capçalera ha de ser la mateixa que el fasta de nucleòtids, afegint el text “, translated by Biopython”.

El fitxer d'entrada que s'ha provat és:

[ls_orchid.fasta](https://github.com/miquelamorosaldev/m14-uf2-intro-biopython/blob/main/ls_orchid.fasta)

El fitxer de sortida generat és:

[ls_orchid_newprot.fasta](https://github.com/miquelamorosaldev/m14-uf2-intro-biopython/blob/main/ls_orchid_newprot.fasta)

## Exercicis bloc 3. Fitxers FASTQ.

Un altre format de fitxer de seqüència molt comú és "FASTQ", que s'utilitza habitualment per informar de dades en experiments de seqüenciació d'alt rendiment:

* Cada entrada del format de fitxer FASTQ comença amb el símbol "@", seguit d'un identificador de seqüència únic. Sovint, aquest identificador codifica informació sobre la màquina i l'execució de seqüenciació a partir de la qual es van produir aquestes dades.
* A continuació està la seqüència de lectura corresponent i després només el símbol "+" o el símbol "+" seguit del mateix identificador únic.
* Finalment, hi ha una cadena de qualitat. La cadena de qualitat té exactament la mateixa longitud que la seqüència de lectura, i cada caràcter codifica una puntuació de qualitat.

Un fitxer d'entrada fastq:

[sequences.fq](https://github.com/miquelamorosaldev/m14-uf2-intro-biopython/blob/main/sequences.fq)

Aquest codi llegeix un fitxer fastq passat per paràmetre des del terminal i genera un gràfic amb la qualitat aproximada: 

[fastqPrintAverageQualityScores.py](https://github.com/miquelamorosaldev/m14-uf2-intro-biopython/blob/main/fastqPrintAverageQualityScores.py)


