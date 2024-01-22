import gzip, os, shutil, urllib3

file ="gene_orthologs"

#####

dir = "data/gene"
if not os.path.isdir(dir):
    os.makedirs(dir)

file_gz = "{}/{}.gz".format(dir,file)
file_tsv = "{}/{}.tsv".format(dir,file)

url = "https://ftp.ncbi.nih.gov/gene/DATA/{}.gz".format(file)

if not os.path.exists(file_tsv):
  
  print("Downloading file {} ...".format(file_gz))
  with urllib3.request("GET", url, preload_content= False) as response:
    
    with open(file_gz, "wb") as writer:
      for chunk in response.stream(10*1024):
        writer.write(chunk)

  print("Extracting file ...")
  with gzip.open(file_gz,"rb") as f_in, open(file_tsv,"wb") as f_out:
      shutil.copyfileobj(f_in, f_out)
  os.remove(file_gz)
