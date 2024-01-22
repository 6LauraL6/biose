from Bio import Entrez
import json
from time import time
from tqdm import tqdm
from urllib.error import HTTPError
from tabulate import tabulate

term = "human[orgn] AND 20[chr] AND alive[prop]"
file = "data/snp_table.txt"

batch_size = 2
Entrez.email = "david@xtec.dev"

with Entrez.esearch("gene", term, retmode="json") as response:
    search = json.loads(response.read())["esearchresult"]


id = search["idlist"][0]
with Entrez.elink(db="snp", dbfrom="gene", linkname="gene_snp", id=id, retmode="json") as response:
    links = json.loads(response.read())
    print(json.dumps(links, indent=2))

'''
The rsID number is a unique label ("rs" followed by a number) used by researchers and databases to identify a specific SNP (Single Nucleotide Polymorphism). It stands for Reference SNP cluster ID and is the naming convention used for most SNPs
'''

'''
count = int(search["count"])
with open(file, "w") as writer:
    for start in tqdm(range(0, count, batch_size), desc=f"Fetch {count} (BS={batch_size}) ..."):
        end = min(count, start+batch_size)

        attempts = 0
        while attempts < 3:
            try:
                with Entrez.elink(db="snp", fromdb="gene", retstart=start, retmax=batch_size, webenv=search["webenv"], query_key=search["querykey"]) as response:
                    links = Entrez.read(response)[0]

                    links = [[link["DbTo"], link["LinkName"], len(
                        link["Link"])] for link in links["LinkSetDb"]]

                    print(tabulate(links, headers=[
                          "DbTo", "LinkName", "len(Link)"]))
                    # links = json.loads(response.read())
                    # print(json.dumps(links,indent=2))
                    quit()
                    # writer.write(data)
                    # break
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    attempts += 1
                    time.sleep(15)
                else:
                    raise

'''
'''
Application 4: Finding unique sets of linked records for each member of a large dataset

Goal: Download separately the SNP rs numbers (identifiers) for each current gene on human chromosome 20.

Solution: First use ESearch to retrieve the Gene IDs for the genes, and then assemble an ELink URL where each Gene ID is submitted as a separate &id parameter.

Input: $query – human[orgn]+AND+20[chr]+AND+alive[prop]

Output: A file named "snp_table" containing on each line the gene id followed by a colon (":") followed by a comma-delimited list of the linked SNP rs numbers.

use LWP::Simple;
use LWP::UserAgent;
$query = 'human[orgn]+AND+20[chr]+AND+alive[prop]';
$db1 = 'gene';
$db2 = 'snp';
$linkname = 'gene_snp';

#assemble the esearch URL
$base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "esearch.fcgi?db=$db1&term=$query&usehistory=y&retmax=5000";

#post the esearch URL
$output = get($url);

#parse IDs retrieved
while ($output =~ /<Id>(\d+?)<\/Id>/sg) {
   push(@ids, $1);
}

#assemble  the elink URL as an HTTP POST call
$url = $base . "elink.fcgi";

$url_params = "dbfrom=$db1&db=$db2&linkname=$linkname";
foreach $id (@ids) {      
   $url_params .= "&id=$id";
}

#create HTTP user agent
$ua = new LWP::UserAgent;
$ua->agent("elink/1.0 " . $ua->agent);

#create HTTP request object
$req = new HTTP::Request POST => "$url";
$req->content_type('application/x-www-form-urlencoded');
$req->content("$url_params");

#post the HTTP request
$response = $ua->request($req); 
$output = $response->content;

open (OUT, ">snp_table") || die "Can't open file!\n";

while ($output =~ /<LinkSet>(.*?)<\/LinkSet>/sg) {

   $linkset = $1;
   if ($linkset =~ /<IdList>(.*?)<\/IdList>/sg) {
      $input = $1;
      $input_id = $1 if ($input =~ /<Id>(\d+)<\/Id>/sg); 
   }

   while ($linkset =~ /<Link>(.*?)<\/Link>/sg) {
      $link = $1;
      push (@output, $1) if ($link =~ /<Id>(\d+)<\/Id>/);
   }
      
   print OUT "$input_id:" . join(',', @output) . "\n";
  
}

close OUT;

Notes: This example uses an HTTP POST request for the elink call, as the number of Gene IDs is over 500. The &retmax parameter in the ESearch call is set to 5000, as this is a reasonable limit to the number of IDs to send to ELink in one request (if you send 5000 IDs, you are effectively performing 5000 ELink operations). If you need to link more than 5000 records, add &retstart to the ESearch call and repeat the entire procedure for each batch of 5000 IDs, incrementing &retstart for each batch.



'''
