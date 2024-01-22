import re


# Utilitzem with open per assegurarnos que es tanca correctament el fitxer.
with open('sequences.fa', 'r') as f:
    lines = f.readlines()


# Aquesta expressió regular agafa les capçaleres
# i les guarda en una llista
hre = re.compile('>(\S+)')


# Aquesta expressió regular agafa les cadenes d'ADN
# i les guarda en una llista.
lre = re.compile('^(\S+)$')


gene = {}
gene_id = None  # Inicialitza gene_id fora del bucle


# Itera sobre les línees del arxiu
for line in lines:
    outh = hre.search(line)
    if outh:
        # Si encuentra un encabezado, guarda el identificador
        gene_id = outh.group(1)
    elif gene_id is not None:
        # Si encuentra una línea de secuencia y tiene un identificador, agrégala al diccionario
        outl = lre.search(line)
        if outl:
            gene[gene_id] = gene.get(gene_id, '') + outl.group(1)


# Imprime los resultados para verificar
for gene_id, sequence in gene.items():
    print(f'ID: {gene_id}, Sequence: {sequence}, Seq.Long: {len(sequence)})')