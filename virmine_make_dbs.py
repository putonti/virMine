import ftputil
import os
import tarfile
import glob
import shutil
import gzip
import csv
from Bio import SeqIO

## VIRAL ##
#retrieve the file from NCBI FTP
print('Getting viral data from NCBI')
host = ftputil.FTPHost('ftp.ncbi.nlm.nih.gov', 'anonymous', 'password')
host.chdir('/genomes/Viruses/')
host.download("all.faa.tar.gz", os.getcwd()+'/all.faa.tar.gz')

output_dir=os.getcwd()+'/v_temp/'
outfilename=os.getcwd()+'/all_viral_faa.fasta'

print('Preparing viral files for database')
#uncompress the file
tar = tarfile.open("all.faa.tar.gz")
for member in tar.getmembers():
  if member.isreg():
    member.name = os.path.basename(member.name)
    tar.extract(member,output_dir)

#concatenate into one file
outfilename=os.getcwd()+'/all_viral_faa.fasta'
with open(outfilename, 'wb') as outfile:
    for filename in glob.glob(output_dir+'*.faa'):
        if filename == outfilename:
            continue
        with open(filename, 'rb') as readfile:
            shutil.copyfileobj(readfile, outfile)

shutil.rmtree(output_dir, ignore_errors=True)

## BACTERIAL ##
print('Getting COG files')
output_dir=os.getcwd()+'/bac_temp/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
host=ftputil.FTPHost('ftp.ncbi.nih.gov','anonymous','password')
host.chdir('/pub/COG/COG2014/data/')
host.download("cog2003-2014.csv",output_dir+'cog2003-2014.csv')
host.download("cognames2003-2014.tab",output_dir+'cognames2003-2014.tab')
host.download("prot2003-2014.fa.gz",output_dir+'prot2003-2014.fa.gz')
with open(output_dir+'/prot2003-2014.fa', 'wb') as f_out, gzip.open(output_dir+'prot2003-2014.fa.gz', 'rb') as f_in:
    shutil.copyfileobj(f_in, f_out)

print('Preparing bacterial files for database. Removing phage/plasmid sequences.')
c=[]
csv_headers=['domain-id','genome-name','protein-id','protein-length','domain-start','domain-end','COG-id','membership-class']
with open(output_dir+'cog2003-2014.csv') as csvfile:
    csv_rows=csv.DictReader(csvfile,csv_headers,delimiter=',')
    for row in csv_rows:
        c.append(row)

t=[]
tab_headers=['COG-id','func','name']
with open(output_dir+'cognames2003-2014.tab') as csvfile:
    tab_rows=csv.DictReader(csvfile,tab_headers,delimiter='\t')
    next(tab_rows, None)
    for row in tab_rows:
        t.append(row)
xs=[l for l in t if l['func']=='X'] #all of the cogs that are phage/plasmid
xCOGids=[i['COG-id'] for i in xs]

xProteinids=[i['protein-id'] for i in c if i['COG-id'] in xCOGids]

output_file=os.getcwd()+'/no_phage_bactTEST.fasta'
output_bact=open(output_file,'wb')
temp_sequences=list(SeqIO.parse(output_dir+'/prot2003-2014.fa','fasta'))
for i in temp_sequences:
    ID=str(i.id)
    proteinID=ID[3:ID.find("|",3)]
    if proteinID in xProteinids:
        continue
    else:
        output_bact.write('>'+str(i.id)+'\n'+str(i.seq)+'\n')
output_bact.close()

shutil.rmtree(output_dir, ignore_errors=True)
