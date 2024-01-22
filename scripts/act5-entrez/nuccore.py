from Bio import Entrez, SeqIO
import os

id = "186972394"
email = "david@xtec.dev"

dir = "data/nuccore"
if not os.path.isdir(dir):
    os.makedirs(dir)

filename = "{}/gi_{}.gbk".format(dir, id)
if not os.path.isfile(filename):
    print("Downloading file ...")
    with Entrez.efetch("nuccore", id=id, rettype="gb", retmode="text", email=email) as response:

        with open(filename, "w") as out_handle:
            out_handle.write(response.read())
        print("Saved file.")

record = SeqIO.read(filename, "genbank")
print(record.seq)
