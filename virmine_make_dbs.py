import ftputil
import os
import tarfile
import glob
import shutil
import gzip
import csv
from Bio import SeqIO

## VIRAL ##
##retrieve the file from NCBI FTP
print('Getting viral data from NCBI')
host = ftputil.FTPHost('ftp.ncbi.nlm.nih.gov', 'anonymous', 'password')
host.chdir('/genomes/Viruses/')
host.download("all.faa.tar.gz", os.getcwd()+'/all.faa.tar.gz')

output_dir=os.getcwd()+'/v_temp/'
outfilename=os.getcwd()+'/all_viral_faa.fasta'

print('Preparing viral files for database')
##uncompress the file
tar = tarfile.open("all.faa.tar.gz")
for member in tar.getmembers():
  if member.isreg():
    member.name = os.path.basename(member.name)
    tar.extract(member,output_dir)

##concatenate into one file
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
output_dir_fix_spaces = output_dir.replace(" ","\ ")

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
host=ftputil.FTPHost('ftp.ncbi.nih.gov','anonymous','password')
host.chdir('/pub/COG/COG2020/data/')
host.download("cog-20.cog.csv",output_dir+'cog-20.cog.csv')
host.download("cog-20.def.tab",output_dir+'cog-20.def.tab')

##OS Command for wget and gunzip
os.system("wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.fa.gz -P "+output_dir_fix_spaces)
os.system("gunzip "+output_dir_fix_spaces+"cog-20.fa.gz")

print('Preparing bacterial files for database. Removing phage/plasmid sequences.')
c=[]
csv_headers=['gene-id','assembly-id','protein-id','protein-length','COG-footprint-coordinates','length-COG-footprint','COG-id','reserved','COG-membership-class', 'PSI-BLAST-bitscore', 'PSI-BLAST-evalue', 'COG-profile-length', 'protein-footprint-coordinates']
with open(output_dir+'cog-20.cog.csv') as csvfile:
    csv_rows=csv.DictReader(csvfile,csv_headers,delimiter=',')
    for row in csv_rows:
        c.append(row)

t=[]
tab_headers=['COG-id','func','name','gene-associated','functional-pathway','pubmed-id','pbd-id']
with open(output_dir+'cog-20.def.tab', mode="r", encoding="cp1252") as csvfile:
    tab_rows=csv.DictReader(csvfile,tab_headers,delimiter='\t')
    for row in tab_rows:
        t.append(row)
xs=[l for l in t if l['func']=='X'] #all of the cogs that are phage/plasmid
xCOGids=[i['COG-id'] for i in xs]

xProteinids=[i['protein-id'] for i in c if i['COG-id'] in xCOGids]

output_file=os.getcwd()+'/no_phage_bact.fasta'
output_bact=open(output_file,'w')
temp_sequences=list(SeqIO.parse(output_dir+'/cog-20.fa','fasta'))
for i in temp_sequences:
    ID=str(i.id)
    if ID[-2] == "_":
        ID = ID[:-2] + '.' + ID[-1:]
    proteinID=ID
    if proteinID in xProteinids:
        continue
    else:
        output_bact.write('>'+str(ID)+'\n'+str(i.seq)+'\n')
output_bact.close()

shutil.rmtree(output_dir, ignore_errors=True)
