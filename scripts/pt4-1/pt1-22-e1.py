from Bio import Entrez
import json
import EntrezUtils
from urllib.error import HTTPError

list_of_terms = ['Homo sapiens tumor protein p53 (TP53), transcript variant 2, mRNA',
                 'Meriones unguiculatus TP53 mRNA for p53, complete cds',
                 'Rattus norvegicus tumor protein p53 (Tp53), mRNA'
#                 'Danio rerio tumor protein p53 (tp53), transcript variant 4, mRNA',
#                 'Gallus gallus tumor protein p53 (TP53), mRNA',
#                 'Canis lupus familiaris p53 gene for P53, complete cds',
#                 'Macaca fascicularis p53 gene, complete cds'
                 ]

for term_el in list_of_terms:
    EntrezUtils.request_search(db='nucleotide',term=term_el,retmax=1,xml_filename=f'./d1/{term_el}.xml')

    species_order = ['homo_sapiens','meriones_unguiculatus','rattus_norvegicus']
#,'rattus_norvegicus','danio_rerio','gallus_gallus','canis_lupus','macaca_fascicularis'
    
    for a, i in zip(list_of_terms, species_order):
        id = EntrezUtils.read_xml(f'./d1/{a}.xml')['IdList'][0]
        name = i
        gb_file = EntrezUtils.request_fetch('nucleotide', id ,'gb',f'./d1/{name}.gb')
        print(gb_file)
        
        