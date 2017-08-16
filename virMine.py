# notes: bbmap, megahit commented out
# for synthetic data, commented out run of sickle b/c they are fasta files


#virMine
import os
import sys
import shutil
import csv
import re
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import collections
import multiprocessing
import argparse
import time
import datetime

##################### FUNCTIONALITY TO GET PARAMETERS FROM COMMAND LINE #####################
def msg(name=None):
    return '''virMine.py [assembly options] [filter options] [databases] -o output_path'''
parser=argparse.ArgumentParser(usage=msg())
groupA=parser.add_mutually_exclusive_group()
groupA.add_argument('-a','--assembler',choices=['spades','metaspades','megahit','all3'],action="store",help='Assembly method. all3 will test all three and pick the best based on N50 scores (Paired-end reads only).')
groupA.add_argument('-A','--assembled_contigs',action="store",metavar='<filename>',help='Contigs file. Start analysis from an existing assembly.')
groupR=parser.add_mutually_exclusive_group()
groupR.add_argument('-s','--single_reads',action="store",metavar='<filename>',help='Single reads.')
groupR.add_argument('-p','--paired_end_reads', nargs=2, metavar=('<filename>','<filename>'), help='Paired-end reads. List both read files.')
parser.add_argument('-o','--output_path',action="store",metavar='<directory>',help='Directory to store all the resulting files (required)')
parser.add_argument('-t','--num_threads',action="store",metavar='<int>',type=int,help='Number of processors to use.')
groupF = parser.add_argument_group('filter options')
groupF.add_argument('-m','--min_contig_size',action="store",metavar='<int>',type=int,help='Minimum contig size.')
groupF.add_argument('-M','--max_contig_size',action="store",metavar='<int>',type=int,help='Maximum contig size.')
groupF.add_argument('-c','--min_coverage',action="store",metavar='<int>',type=int,help='Minimum coverage.')
groupF.add_argument('-cov','--min_SPAdes_cov',action="store",metavar='<float>',type=float,help='Minimum SPAdes cov value.')
groupF.add_argument('-g','--genes_of_interest',action="store",metavar='<filename>',help='PROTEIN sequences of interest. FASTA format. Contig must contain homology to genes of interest.')
groupDB = parser.add_argument_group('database options')
groupDB.add_argument('-v','--viral',action="store",metavar='<filename>',help='Viral PROTEIN sequence collection. FASTA format. (required)')
groupDB.add_argument('-nv','--non_viral',action="store",metavar='<filename>',help='Non-viral PROTEIN sequence collection. FASTA format. (required)')
parser.add_argument('-b','--blast_option',choices=['blastx','blastp'],action="store",help='BLAST method to classify as viral or non-viral. Note: blastp is faster.')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args=parser.parse_args()
if args.assembler is None and args.assembled_contigs is None:
    parser.error('Assembler method or assembly file must be specified.')
if args.assembler is not None:
    if args.paired_end_reads is None and args.single_reads is None:
        parser.error('Reads must be provided for assembly.')
if args.output_path is None:
    parser.error('Output path must be provided.')
if args.assembler=='all3' and args.single_reads is not None:
    parser.error('The all3 method can only examine paired-end reads.')
if args.assembler=='all3' and args.min_SPAdes_cov is not None:
    print('Cannot guarantee SPAdes or metaSPAdes will be the best assembly. Ignoring cov filter.')
    args.min_SPAdes_cov=None
if args.assembler=='metaspades' and args.single_reads is not None:
    parser.error('metaSPAdes cannot process single end reads.')
if args.assembler=='megahit' and args.min_SPAdes_cov is not None:
    print('Megahit assembly can not be filtered by SPAdes cov values. Ignoring cov filter.')
    args.min_SPAdes_cov=None
if args.viral is None or args.non_viral is None:
    parser.error('You must specify both the viral and nonviral databases.')
##################### FUNCTIONALITY TO GET PARAMETERS FROM COMMAND LINE #####################

#PROGRAM FUNCTIONALITY
#CHECK TO MAKE CERTAIN ALL OF THE PATHS ARE SET
def path_check(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

#CHECK ALL PATHS
def check_installs():
    #check sickle
    if path_check('sickle') == None:
        print('Program SICKLE not in path')
        return False

    #check SPAdes
    if path_check('spades.py') == None:
        print('Program SPAdes not in path')
        return False

    #check MEGAHIT
#    if path_check('megahit') == None:
#        print('Program MEGAHIT not in path')
#        return False

    #check GLIMMER
#    if path_check('glimmer3') == None or path_check('long-orfs') == None or path_check('extract') == None or path_check('build-icm') == None:
#        print ('Program GLIMMER not in path')
#        return False

    #check BLASTX
    if path_check('blastx') == None:
        print('BLAST+ tools are not added to the path')

    #check BBMAP
#    if path_check('bbmap.sh') == None:
#        print('Program bbmap not in path')
#        return False

   
#CHECK IF A FILE EXISTS
def check_for_file(file_name):
    if not os.path.isfile(file_name):
        return False
    else:
        return True

#CHECK IF A GROUP OF FILES EXISTS
def check_files(*file_names):
    error_message=''
    for i in file_names:
        if check_for_file(i)==False:
            error_message+='Error: File '+i+' does not exist.\n'
    if len(error_message)==0:
        return True,error_message
    else:
        return False,error_message

#write to log
def log(message,file_name):
    logf=open(file_name,'a')
    logf.write(message+'\n')
    logf.close()

#GET NUMBER OF SEQUENCES IN A FASTA FORMAT FILE
def get_number_of_sequences(contig_file):
    temp_sequences=list(SeqIO.parse(contig_file,'fasta'))
    num_seqs=len(temp_sequences)
    del temp_sequences[:]
    return num_seqs

#GENERATE A LIST OF THE CONTIG HEADERS
def get_list_of_contigs(contig_file):
    contigs=[]
    filename=open(contig_file,'r')
    for l in filename:
        if l[0]=='>':
            contigs.append(l[1:].strip())
    filename.close()
    return contigs

#PARSE THE BLAST AND RETURN A LIST OF DICTIONARIES
def parse_blast(filename,headers):
    x=[]
    blast_results=open(filename,'r')
    rows=csv.DictReader(blast_results,headers,delimiter=',')
    for row in rows:
        x.append(row)
    blast_results.close()
    return x

#WRITE OUT BLAST RESULTS (csv file) FOR CONTIGS OF ITNEREST
def write_out_blast_results(blast_input,output_file,contigs_of_interest):
    headers=['qseqid','sseqid','qcovs','qstart','qend','pident','length','evalue','bitscore']
    x=parse_blast(blast_input,headers)
    y=[]
    for i in x:
        contig=i['qseqid'][:i['qseqid'].find('_orf')]
        if contig in contigs_of_interest:
                y.append(i)

    outputf=open(output_file,'w')
    outputf.write('contig orf,hit,start (query),stop (query),query coverage,%ID,length,bitscore\n')
    for i in y:
        outputf.write(i['qseqid']+','+i['sseqid']+','+i['qstart']+','+i['qend']+','+i['qcovs']+','+i['pident']+','+i['length']+','+i['bitscore']+'\n')
    outputf.close()
    return True

#WRITE OUT FASTA FORMAT SEQUENCES FOR CONTIGS OF INTEREST
def write_out_contigs_of_interest(contig_file,output_file,contigs_of_interest):
    sequences=list(SeqIO.parse(contig_file,'fasta'))
    outp=open(output_file,'w')
    for i in sequences:
        if i.description in contigs_of_interest:
            outp.write('>'+i.description+'\n')
            outp.write(str(i.seq)+'\n')
    outp.close()

#WRITE OUT FASTA FORMAT ORF SEQUENCES FOR CONTIGS OF INTEREST
def write_out_orfs_for_contigs_of_interest(orf_file,output_file,contigs_of_interest):
    sequences=list(SeqIO.parse(orf_file,'fasta'))
    outp=open(output_file,'w')
    for i in sequences:
        contig_name=i.description[:i.description.find('_orf')]
        if contig_name in contigs_of_interest:
            outp.write('>'+i.description+'\n')
            outp.write(str(i.seq)+'\n')
    outp.close()

#RUN BLAST AGAINST VIRAL AND NON-VIRAL DATABASES --> blastx
#blast input must have the contig and orf information within the header. no spaces.... note orf files must be preprocessed (DONE IN RUN_GLIMMER)
def run_blasts(orf_file,bact_db_file,viral_db_file,output_path,num_procs):
    print('Identifying non-viral coding regions within contigs')
    headers=['qseqid','sseqid','qcovs','qstart','qend','pident','length','evalue','bitscore']
    bact_blastx_cline=NcbiblastxCommandline(query=orf_file, db=bact_db, num_threads=num_procs, max_target_seqs=1, outfmt='"10 qseqid sseqid qcovs qstart qend pident length evalue bitscore"', out=output_path+'/all_nonviral.blastx')
    stdout, stderr = bact_blastx_cline()
    b=parse_blast(output_path+'/all_nonviral.blastx',headers)

    print('Identifying viral coding regions within contigs')
    viral_blastx_cline=NcbiblastxCommandline(query=orf_file, db=viral_db, num_threads=num_procs, max_target_seqs=1, outfmt='"10 qseqid sseqid qcovs qstart qend pident length evalue bitscore"', out=output_path+'/all_viral.blastx')
    stdout, stderr = viral_blastx_cline()
    v=parse_blast(output_path+'/all_viral.blastx',headers)
    return b,v

def make_protein_record(nuc_record):
    return SeqRecord(seq=nuc_record.seq.translate(to_stop=True),id=nuc_record.id,description=nuc_record.description)

#RUN BLAST AGAINST VIRAL AND NON-VIRAL DATABASES --> blastp
#blast input must have the contig and orf information within the header. no spaces.... note orf files must be preprocessed (DONE IN RUN_GLIMMER)
def run_blasts_p(orf_file,bact_db_file,viral_db_file,output_path,num_procs):
    #make aa file for orfs
    proteins=(make_protein_record(nuc_rec) for nuc_rec in SeqIO.parse(orf_file,'fasta'))
    aa_orf_file=output_path+'/temp/aa_predicted_orfs.fasta'
    SeqIO.write(proteins,aa_orf_file,'fasta')
    
    print('Identifying non-viral coding regions within contigs')
    headers=['qseqid','sseqid','qcovs','qstart','qend','pident','length','evalue','bitscore']
    bact_blastp_cline=NcbiblastpCommandline(query=aa_orf_file, db=bact_db, num_threads=num_procs, max_target_seqs=1, outfmt='"10 qseqid sseqid qcovs qstart qend pident length evalue bitscore"', out=output_path+'/all_nonviral.blastp')
    stdout, stderr = bact_blastp_cline()
    b=parse_blast(output_path+'/all_nonviral.blastp',headers)

    print('Identifying viral coding regions within contigs')
    viral_blastp_cline=NcbiblastpCommandline(query=aa_orf_file, db=viral_db, num_threads=num_procs, max_target_seqs=1, outfmt='"10 qseqid sseqid qcovs qstart qend pident length evalue bitscore"', out=output_path+'/all_viral.blastp')
    stdout, stderr = viral_blastp_cline()
    v=parse_blast(output_path+'/all_viral.blastp',headers)
    return b,v

#DETERMINE IF CONTIG IS VIRAL OR NON-VIRAL --> blastx
def determine_contig_scores(contig_file,orf_file,bact_db,viral_db,output_path,num_threads):
    b,v=run_blasts(orf_file,bact_db,viral_db,output_path,num_threads)
    viral_scores=collections.defaultdict(int)
    bact_scores=collections.defaultdict(int)
    contigs=get_list_of_contigs(contig_file)
    for i in contigs:
        viral_scores[i]=0
        bact_scores[i]=0

    #contig name is only part of the orf_file header
    #assumes orf_file header is >contigname_orfID_start_stop_len
    for i in v:
        contig=i['qseqid'][:i['qseqid'].find('_orf')]
        found_match=False
        for j in b:
            if i['qseqid']==j['qseqid']:
                if float(i['bitscore'])>float(j['bitscore']):
                    viral_scores[contig]+=float(i['bitscore'])
                    found_match=True
                else:
                    bact_scores[contig]+=float(j['bitscore'])
                    found_match=True
        if found_match==False:
            viral_scores[contig]+=float(i['bitscore'])
            
    for i in b:
        contig=i['qseqid'][:i['qseqid'].find('_orf')]
        found_match=False
        for j in v:
            if i['qseqid']==j['qseqid']:
                found_match=True
        if found_match==False:
            bact_scores[contig]+=float(i['bitscore'])

    num_called_bacteria=0
    num_called_viral=0
    num_called_unkn=0
    v_contigs=[]
    unkn_contigs=[]
    for i in contigs:
        if viral_scores[i]>bact_scores[i]:
            v_contigs.append(i)
            num_called_viral+=1
        else:
            if viral_scores[i]==0 and bact_scores[i]==0:
                unkn_contigs.append(i)
                num_called_unkn+=1
            else:
                num_called_bacteria+=1
    return num_called_bacteria, num_called_viral, num_called_unkn, v_contigs, unkn_contigs

#DETERMINE IF CONTIG IS VIRAL OR NON-VIRAL --> blastp
def determine_contig_scores_p(contig_file,orf_file,bact_db,viral_db,output_path,num_threads):
    b,v=run_blasts_p(orf_file,bact_db,viral_db,output_path,num_threads)
    viral_scores=collections.defaultdict(int)
    bact_scores=collections.defaultdict(int)
    contigs=get_list_of_contigs(contig_file)
    for i in contigs:
        viral_scores[i]=0
        bact_scores[i]=0

    #contig name is only part of the orf_file header
    #assumes orf_file header is >contigname_orfID_start_stop_len
    for i in v:
        contig=i['qseqid'][:i['qseqid'].find('_orf')]
        found_match=False
        for j in b:
            if i['qseqid']==j['qseqid']:
                if float(i['bitscore'])>float(j['bitscore']):
                    viral_scores[contig]+=float(i['bitscore'])
                    found_match=True
                else:
                    bact_scores[contig]+=float(j['bitscore'])
                    found_match=True
        if found_match==False:
            viral_scores[contig]+=float(i['bitscore'])
            
    for i in b:
        contig=i['qseqid'][:i['qseqid'].find('_orf')]
        found_match=False
        for j in v:
            if i['qseqid']==j['qseqid']:
                found_match=True
        if found_match==False:
            bact_scores[contig]+=float(i['bitscore'])

    num_called_bacteria=0
    num_called_viral=0
    num_called_unkn=0
    v_contigs=[]
    unkn_contigs=[]
    for i in contigs:
        if viral_scores[i]>bact_scores[i]:
            v_contigs.append(i)
            num_called_viral+=1
        else:
            if viral_scores[i]==0 and bact_scores[i]==0:
                unkn_contigs.append(i)
                num_called_unkn+=1
            else:
                num_called_bacteria+=1
    return num_called_bacteria, num_called_viral, num_called_unkn, v_contigs, unkn_contigs

#REFORMATTING OF THE CONTIG HEADERS
def qc_contig_file(contig_file):
    sequences=list(SeqIO.parse(contig_file,'fasta'))
    outputfile=open(contig_file,'w')
    for i in sequences:
        description=i.description
        description=re.sub(r'\,','',description)
        description=re.sub(r'\ ','_',description)
        outputfile.write('>'+description+'\n')
        outputfile.write(str(i.seq)+'\n')
    outputfile.close()
    return contig_file

#RUN SICKLE
def run_sickle(output_path,r_type='pe',*read_files):
    sickle_command='sickle '
    if r_type=='pe':    #paired-end reads
        if len(read_files)!=2:
            print('Paired-end reads requires 2 read files entered.')
            return False
        else:
            reads_f=read_files[0]
            reads_r=read_files[1]
            output_file_label1='trimmed_'+read_files[0][read_files[0].rfind("/")+1:]
            output_file_label2='trimmed_'+read_files[1][read_files[0].rfind("/")+1:]
            sickle_command+='pe -f '+read_files[0]+' -r '+read_files[1]+' -t sanger -o '+output_path+'/temp/'+output_file_label1\
                         +' -p '+output_path+'/temp/'+output_file_label2+' -s '+output_path+'/temp/singletons.fastq -l 30'
            check_output=output_path+'/temp/'+output_file_label1
            print(sickle_command)
    else:
        if len(read_files)!=1:
            print('Single-end reads requires 1 read file entered.')
            return False
        else:
            read1=read_files[0]
            output_file_label='trimmed_'+read1[read1.rfind("/")+1:]
            sickle_command+='se -f '+read1+' -t sanger -o '+output_path+'/temp/'+output_file_label+' -l 100'
            check_output=output_path+'/temp/'+output_file_label
            print(sickle_command)
    os.system(sickle_command)
    if os.path.getsize(check_output)>0:
        return True
    else:
        return False
    
#RUN SPAdes
def run_spades(num_threads,output_path,r_type='pe',*read_files):
    spades_command='spades.py'
    if r_type=='pe':

        if len(read_files)!=2:
            print('Paired-end reads requires 2 read files entered.')
            return False
        else:
            spades_command+=' -k 33,55,77,99,127 -t '+num_threads+' --only-assembler -1 '+read_files[0]+' -2 '+read_files[1]+' -o '+output_path+'/temp/assembly/'
            print(spades_command)
    else:
        if len(read_files)!=1:
            print('Single-end reads requires 1 read file entered.')
            return False
        else:
            spades_command+=' -k 33,55,77,99,127 -t '+num_threads+' --only-assembler -s '+read_files[0]+' -o '+output_path+'/temp/assembly/'
            print(spades_command)
    os.system(spades_command)
    return check_for_file(output_path+'/temp/assembly/scaffolds.fasta')

#RUN metaSPAdes
def run_metaspades(num_threads,output_path,r_type='pe',*read_files):
    metaspades_command='spades.py'
    if r_type=='pe':
        if len(read_files)!=2:
            print('Paired-end reads requires 2 read files entered.')
            log('Paired-end reads requires 2 read files entered.',output_path+'/virMine_log.txt')
            return False
        else:
            metaspades_command+=' -k 33,55,77,99,127 --meta -t '+num_threads+' --only-assembler -1 '+read_files[0]+' -2 '+read_files[1]+' -o '+output_path+'/temp/assembly/'
            print(metaspades_command)
    else:
        print('metaSPAdes does not presently work with single end reads.')
        log('metaSPAdes does not presently work with single end reads.',output_path+'/virMine_log.txt')
        return False
    os.system(metaspades_command)
    return check_for_file(output_path+'/temp/assembly/scaffolds.fasta')

#RUN MEGAHIT
def run_megahit(num_threads,output_path,r_type='pe',*read_files):
    megahit_command='megahit -t '+num_threads
    if r_type=='pe':
        if len(read_files)!=2:
            print('Paired-end reads requires 2 read files entered.')
            return False
        else:
            megahit_command+=' -1 '+read_files[0]+' -2 '+read_files[1]+' -o '+output_path+'/temp/assembly/'
            print(megahit_command)
    else:
        if len(read_files)!=1:
            print('Single-end reads requires 1 read file entered.')
            return False
        else:
            megahit_command+=' -r '+read_files[0]+' -o '+output_path+'/temp/assembly/'
            print(megahit_command)
    os.system(megahit_command)
    return True

#RUN ALL ASSEMBLY OPTIONS
def run_all_assemblers(num_threads,output_path,r_type='pe',*read_files):
    log('Running spades, metaspades & megahit assemblies.',output_path+'/virMine_log.txt')##
    spades_command='spades.py -k 33,55,77,99,127 -t '+num_threads+' --only-assembler -1 '+read_files[0]+' -2 '+read_files[1]+' -o '+output_path+'/temp/spades/assembly/'
    os.system(spades_command)
    metaspades_command='spades.py -k 33,55,77,99,127 --meta -t '+num_threads+' --only-assembler -1 '+read_files[0]+' -2 '+read_files[1]+' -o '+output_path+'/temp/metaspades/assembly/'
    os.system(metaspades_command)
    os.makedirs(output_path+'/temp/megahit') #make megahit directory
    megahit_command='megahit -t '+num_threads+' -1 '+read_files[0]+' -2 '+read_files[1]+' -o '+output_path+'/temp/megahit/assembly/'
    print(megahit_command)
    os.system(megahit_command)
    return True

#CALCULATE N50
def N50(file_name):
    sequences=list(SeqIO.parse(file_name,'fasta'))
    list_of_lens=[]
    sum_of_all=0
    for i in sequences:
        list_of_lens.append(len(i.seq))
        sum_of_all+=len(i.seq)
    list_of_lens.sort(reverse=True)
    halfsies=float(sum_of_all)/2.0
    val=0
    for i in list_of_lens:
        val+=i
        if val>=halfsies:
            return i

#FIND BEST ASSEMBLY AND PUT IN /TEMP/ASSEMBLY
def find_best(output_path):
    os.makedirs(output_path+'/temp/assembly') #make temp/assembly directory
    if check_for_file(output_path+'/temp/spades/assembly/scaffolds.fasta')==True:
        spades_score=N50(output_path+'/temp/spades/assembly/scaffolds.fasta')
    else:
        spades_score=0
    if check_for_file(output_path+'/temp/metaspades/assembly/scaffolds.fasta')==True:
        metaspades_score=N50(output_path+'/temp/metaspades/assembly/scaffolds.fasta')
    else:
        metaspades_score=0
    if check_for_file(output_path+'/temp/megahit/assembly/final.contigs.fa')==True:
        megahit_score=N50(output_path+'/temp/megahit/assembly/final.contigs.fa')
    else:
        megahit_score=0
    if spades_score==metaspades_score==megahit_score==0:
        print('Problem with assemblies.')
        log('Assemblies did not run.',output_path+'/virMine_log.txt')##
        return False
    best_alignment=''
    if int(spades_score) >= int(metaspades_score) and int(spades_score) >= int(megahit_score):
        best_alignment='spades'
        src=output_path+'/temp/spades/assembly/scaffolds.fasta'
    if int(metaspades_score) >= int(spades_score) and int(metaspades_score) >= int(megahit_score):
        best_alignment='metaspades'
        src=output_path+'/temp/metaspades/assembly/scaffolds.fasta'
    if int(megahit_score) >= int(spades_score) and int(megahit_score) >= int(metaspades_score):
        best_alignment='megahit'
        src=output_path+'/temp/megahit/assembly/final.contigs.fa'
    log('Assembly scores:\n\tspades: '+str(spades_score)+'\n\tmetaspades: '+str(metaspades_score)+'\n\tmegahit: '+str(megahit_score)+'\nBest alignment: '+best_alignment,output_path+'/virMine_log.txt')##
    shutil.copy(src,output_path+'/temp/assembly/scaffolds.fasta') #copy assembly into /temp/assembly
    return True

#RUN GLIMMER
def run_glimmer(contig_file,output_path):
    #predict the genes
    glimmer_command='./g3-iterated-viral.csh '+contig_file+' '+output_path+'/temp/orfs'
    print(glimmer_command)
    os.system(glimmer_command)

    if check_for_file(output_path+'/temp/orfs.coords')==True:
        if os.path.getsize(output_path+'/temp/orfs.coords')==0:
            print('GLIMMER did not predict any genes.')
            log('GLIMMER did not predict any genes.',output_path+'/virMine_log.txt')##
            return False
    else:
        print('GLIMMER did not predict any genes.')
        log('GLIMMER did not predict any genes.',output_path+'/virMine_log.txt')##
        return False        

    #determine which orfs are in which contigs
    contigs=get_list_of_contigs(contig_file)
    contig_label=contigs[0]
    orfs=collections.OrderedDict()
    orfies=open(output_path+'/temp/orfs.coords','r')
    temp_orfs=[i.strip().split() for i in orfies]
    for i in temp_orfs:
        if i[0][0]=='>': #new contig
            contig_label=i[0][1:]
        else:
            orfs[(i[0],contig_label)]=[i[1],i[2],i[3]]
    orfies.close()

    #multi-extract function -- includes the last 3 nucleotides (the stop codon)
    sequences=list(SeqIO.parse(contig_file,'fasta'))
    orf_sequences=[]
    for i,j in orfs: #i is orfID, j is the contig
        for k in sequences:   
            if j==k.description:
                start,stop,frame=orfs[(i,j)]
                start=int(start)
                stop=int(stop)
                frame=int(frame)
                if start>stop:
                    if frame<0: #it's reverse complement
                        subseq=k.seq[stop-1:start].reverse_complement()
                    else: #it goes over the origin
                        subseq=k.seq[start-1:]+k.seq[:stop]
                else:
                    if frame<0: #it goes over the origin and reverse complement
                        subseq=(k.seq[stop-1:]+k.seq[:start]).reverse_complement()
                    else:
                        subseq=k.seq[start-1:stop]
                orf_label='>'+j+'_'+i+'_'+str(start)+'_'+str(stop)
                orf_sequences.append([orf_label,subseq])
                break
    orf_output=open(output_path+'/predicted_orfs.fasta','w')
    for i in orf_sequences:
        orf_output.write(str(i[0])+'\n')
        orf_output.write(str(i[1])+'\n')
    orf_output.close()
    log('ORF sequences written to file: '+output_path+'/predicted_orfs.fasta',output_path+'/virMine_log.txt')##
    return output_path+'/predicted_orfs.fasta'

#RUN READ QC (SICKLE) AND ASSEMBLY (SPADES, METASPADES, OR MEGAHIT)
def qc_and_assembly(assembler,num_threads,type_reads,*reads):
    if assembler=='spades':
        if type_reads=='pe':
            read1=reads[0][0]
            read2=reads[0][1]
            if run_sickle(output_path,'pe',read1,read2)==True:
                log('Read QC completed by sickle.',output_path+'/virMine_log.txt')##
                output_file_label1='trimmed_'+read1[read1.rfind("/")+1:]
                output_file_label2='trimmed_'+read2[read2.rfind("/")+1:]
                if run_spades(num_threads,output_path,'pe',output_path+'/temp/'+output_file_label1,output_path+'/temp/'+output_file_label2)==True:
                    log('SPAdes assembly written to: '+output_path+'/temp/assembly/',output_path+'/virMine_log.txt')##
                    return True
                else:
                    return False
            else:
                print('There was an error with read QC. Check files and try again.\n')
                log('There was an error with read QC. Check files and try again.',output_path+'/virMine_log.txt')##
                return False
        else:
            read1=reads[0]
            if True: #run_sickle(output_path,'se',read1)==True:                                                            ########## remove when not synthetic
                log('Read QC completed by sickle.',output_path+'/virMine_log.txt')##
                output_file_label='trimmed_'+read1[read1.rfind("/")+1:]
                os.system('cp '+read1+' '+output_path+'/temp/'+output_file_label)                                          ########## remove when not synthetic
                if run_spades(num_threads,output_path,'se',output_path+'/temp/'+output_file_label)==True:
                    log('SPAdes assembly written to: '+output_path+'/temp/assembly/',output_path+'/virMine_log.txt')##
                    return True
                else:
                    return False
            else:
                print('There was an error with read QC. Check files and try again.\n')
                log('There was an error with read QC. Check files and try again.',output_path+'/virMine_log.txt')##
                return False
    elif assembler=='metaspades':
        if type_reads=='pe':
            read1=reads[0][0]
            read2=reads[0][1]
            if run_sickle(output_path,'pe',read1,read2)==True:
                log('Read QC completed by sickle.',output_path+'/virMine_log.txt')##
                output_file_label1='trimmed_'+read1[read1.rfind("/")+1:]
                output_file_label2='trimmed_'+read2[read2.rfind("/")+1:]
                if run_metaspades(num_threads,output_path,'pe',output_path+'/temp/'+output_file_label1,output_path+'/temp/'+output_file_label2)==True:
                    log('metaSPAdes assembly written to: '+output_path+'/temp/assembly/',output_path+'/virMine_log.txt')##
                    return True
                else:
                    return False
            else:
                print('There was an error with read QC. Check files and try again.\n')
                log('There was an error with read QC. Check files and try again.',output_path+'/virMine_log.txt')##
                return False
        else:
            read1=reads[0]
            if run_sickle(output_path,'se',read1)==True:
                log('Read QC completed by sickle.',output_path+'/virMine_log.txt')##
                output_file_label='trimmed_'+read1[read1.rfind("/")+1:]
                if run_metaspades(num_threads,output_path,'se',output_path+'/temp/'+output_file_label)==True:
                    log('metaSPAdes assembly written to: '+output_path+'/temp/assembly/',output_path+'/virMine_log.txt')##
                    return True
                else:
                    return False
            else:
                print('There was an error with read QC. Check files and try again.\n')
                log('There was an error with read QC. Check files and try again.',output_path+'/virMine_log.txt')##
                return False
    elif assembler=='megahit':
        if type_reads=='pe':
            read1=reads[0][0]
            read2=reads[0][1]
            if run_sickle(output_path,'pe',read1,read2)==True:
                log('Read QC completed by sickle.',output_path+'/virMine_log.txt')##
                output_file_label1='trimmed_'+read1[read1.rfind("/")+1:]
                output_file_label2='trimmed_'+read2[read2.rfind("/")+1:]
                if run_megahit(num_threads,output_path,'pe',output_path+'/temp/'+output_file_label1,output_path+'/temp/'+output_file_label2)==True:
                    log('megahit assembly written to: '+output_path+'/temp/assembly/',output_path+'/virMine_log.txt')##
                    return True
                else:
                    return False
            else:
                print('There was an error with read QC. Check files and try again.\n')
                log('There was an error with read QC. Check files and try again.',output_path+'/virMine_log.txt')##
                return False
        else:
            read1=reads[0]
            if run_sickle(output_path,'se',read1)==True:
                log('Read QC completed by sickle.',output_path+'/virMine_log.txt')##
                output_file_label='trimmed_'+read1[read1.rfind("/")+1:]
                if run_megahit(num_threads,output_path,'se',output_path+'/temp/'+output_file_label)==True:
                    log('megahit assembly written to: '+output_path+'/temp/assembly/',output_path+'/virMine_log.txt')##
                    return True
                else:
                    return False
            else:
                print('There was an error with read QC. Check files and try again.\n')
                log('There was an error with read QC. Check files and try again.',output_path+'/virMine_log.txt')##
                return False
    elif assembler=='all3':
        if type_reads=='pe':
            read1=reads[0][0]
            read2=reads[0][1]
            if run_sickle(output_path,'pe',read1,read2)==True:
                log('Read QC completed by sickle.',output_path+'/virMine_log.txt')##
                output_file_label1='trimmed_'+read1[read1.rfind("/")+1:]
                output_file_label2='trimmed_'+read2[read2.rfind("/")+1:]
                #run all 3 assemblies and pick the best
                if run_all_assemblers(num_threads,output_path,'pe',output_path+'/temp/'+output_file_label1,output_path+'/temp/'+output_file_label2)==True:
                    return find_best(output_path)
            else:
                print('There was an error with read QC. Check files and try again.\n')
                log('There was an error with read QC. Check files and try again.',output_path+'/virMine_log.txt')##
                return False
        else:
            print('The all3 method requires paired-end reads. Assembly not performed.\n')
            log('The all3 method requires paired-end reads. Assembly not performed.',output_path+'/virMine_log.txt')##
            return False

#FILTER OPTIONS m AND M
# filter contigs by size
def filter_contigs_by_size(contig_file,condition,size):
    counter=0
    sequences=list(SeqIO.parse(contig_file,'fasta'))
    outputfile=open(contig_file,'w')
    for i in sequences:
        if condition=='>':
            if len(i.seq)>=int(size):
                outputfile.write('>'+i.description+'\n')
                outputfile.write(str(i.seq)+'\n')
                counter=counter+1
        else:
            if len(i.seq)<=int(size):
                outputfile.write('>'+i.description+'\n')
                outputfile.write(str(i.seq)+'\n')
                counter=counter+1
    outputfile.close()
    log('Filter by size... '+condition+str(size)+': '+str(counter),output_path+'/virMine_log.txt')##
    del sequences[:]
    return contig_file

#FILTER OPTION C
def filter_by_coverage(contig_file,val,type_reads,*read_files):
    sequences=list(SeqIO.parse(contig_file,'fasta'))
    log('Filter by coverage... >='+str(val)+': '+str(len(sequences)),output_path+'/virMine_log.txt')##
    if type_reads=='pe':       
        read1=read_files[0]
        read2=read_files[1]
        message='Read files to calculate coverage: '+read1+' & '+read2
        log(message,output_path+'/virMine_log.txt')##
        bbwrap_command='bbwrap.sh ref='+contig_file+' in='+read1+' in2='+read2+' out='+output_path+'/temp/aln.sam.gz kfilter=22 subfilter=15 maxindel=80'
        os.system(bbwrap_command)
        print(bbwrap_command)
        if check_for_file(output_path+'/temp/aln.sam.gz')==False:
            print('Error with BBMap. Filter by coverage not performed.')
            log('Error with BBMap. Filter by coverage not performed.',output_path+'/virMine_log.txt')##
            return contig_file
        pileup_command='pileup.sh in='+output_path+'/temp/aln.sam.gz '+'out='+output_path+'/temp/covstats.txt twocolumn=t'
        os.system(pileup_command)
        print(pileup_command)
        if check_for_file(output_path+'/temp/covstats.txt')==False:
            print('Error with pileup. Filter by coverage not performed.')
            log('Error with pileup. Filter by coverage not performed.',output_path+'/virMine_log.txt')##
            return contig_file            
    else:
        message='Read file to calculate coverage: '+read_files[0],output_path+'/virMine_log.txt'
        log(message,output_path+'/virMine_log.txt')##
        bbwrap_command='bbwrap.sh ref='+contig_file+' in='+read_files[0]+' out='+output_path+'/temp/aln.sam.gz kfilter=22 subfilter=15 maxindel=80'
        os.system(bbwrap_command)
        print(bbwrap_command)
        if check_for_file(output_path+'/temp/aln.sam.gz')==False:
            print('Error with BBMap. Filter by coverage not performed.')
            log('Error with BBMap. Filter by coverage not performed.',output_path+'/virMine_log.txt')##
            return contig_file
        pileup_command='pileup.sh in='+output_path+'/temp/aln.sam.gz '+'out='+output_path+'/temp/covstats.txt twocolumn=t'
        os.system(pileup_command)
        print(pileup_command)
        if check_for_file(output_path+'/temp/covstats.txt')==False:
            print('Error with pileup. Filter by coverage not performed.')
            log('Error with pileup. Filter by coverage not performed.',output_path+'/virMine_log.txt')##
            return contig_file    
    #PARSE and write to file only those contigs with a coverage value >= val
    val=float(val)
    meet_threshold=[]
    cov_file=output_path+'/temp/covstats.txt'
    with open(cov_file) as csvfile:
        reader=csv.DictReader(csvfile,delimiter='\t')
        for row in reader:
            if float(row['Avg_fold'])>=val:
                meet_threshold.append(row)
    outputfile=open(contig_file,'w')
    for i in meet_threshold:
        outputfile.write('>'+i['#ID']+'\n')
        for j in sequences:
            if i['#ID']==j.id:
                outputfile.write(str(j.seq)+'\n')
                break
    outputfile.close()
    log('Filter by cov... >='+str(val)+': '+str(len(meet_threshold)),output_path+'/virMine_log.txt')##
    return contig_file

#FILTER OPTION COV
# filter contigs by spades cov value
def filter_by_spades_cov(contig_file,val):
    val=float(val)
    meet_threshold=[]
    no_cov_val=False
    sequences=list(SeqIO.parse(contig_file,'fasta'))
    for i in sequences:
        try:
            covScore=i.id[i.id.index('_cov_'):]
            covScore=covScore[5:]
            covScore=float(covScore)
            if covScore>=val:
                meet_threshold.append(i)
        except:
            no_cov_val=True
    if no_cov_val==True:
        print('cov value not found in contigs. Contigs were not filtered.')
    else:
        outputfile=open(contig_file,'w')
        for i in meet_threshold:
            outputfile.write('>'+i.description+'\n')
            outputfile.write(str(i.seq)+'\n')
        outputfile.close()
        log('Filter by cov... >='+str(val)+': '+str(len(meet_threshold)),output_path+'/virMine_log.txt')##
    del sequences[:]
    del meet_threshold[:]
    return contig_file

#FILTER OPTION G
# keep only those showing homology to sequences of interest
def filter_by_seq_of_interest(contig_file,file_name,num_procs):
    #make database of sequences of interest
    db_genes=output_path+'/temp/genes_of_interest'
    output_file=output_path+'/temp/seqs_of_interest.blastx'
    makeblast_command='makeblastdb -in '+file_name+' -out '+db_genes+' -title '+db_genes+' -dbtype prot'
    os.system(makeblast_command)
    fileQC,errmsg=check_files(db_genes+'.pin',viral_db+'.pin')
    if fileQC==False:
        print('Error creating database of genes of interest.\nFilter skipped.\n')
        log('Error creating database of genes of interest. Filter skipped. '+errmsg,output_path+'/virMine_log.txt')##
        return contig_file
    else:
        log('Filter by genes of interest... :\n\tDatabase of genes of interest: '+db_genes,output_path+'/virMine_log.txt')##
        blastx_cline=NcbiblastxCommandline(query=contig_file, db=db_genes, num_threads=num_procs, max_target_seqs=1, outfmt='"10 qseqid sseqid bitscore"', out=output_file)
        stdout, stderr = blastx_cline()
        headers=['qseqid','sseqid','bitscore']
        b=parse_blast(output_file,headers)
        keep=[]
        for i in b:
            if i['bitscore']>50:
                keep.append(i)
        del b[:]
        if(len(keep)>0):
            sequences=list(SeqIO.parse(contig_file,'fasta'))
            outputfile=open(contig_file,'w')
            for i in keep:
                for j in sequences:
                    if i['qseqid']==j.description:
                        outputfile.write('>'+j.description+'\n')
                        outputfile.write(str(j.seq)+'\n')
            outputfile.close()
            log('\tNumber of hits to genes of interest: '+str(len(keep)),output_path+'/virMine_log.txt')##
            del sequences[:]
            del keep[:]
        else:
            print('No contigs share homology to the sequences of interest. Filter ignored.')
        return contig_file

#FILTERS THE CONTIGS ACCORDING TO THE FILTER OPTIONS
# current options: min_size, max_size, min_coverage, min_cover_spades, genes_of_interest
# filter parameters are listed in variable filter_params separated by pipe
# parameter codes m=min_size; M=max_size; c=coverage; s=min_cover_spades; g=fasta sequence(s) of interest
def filter_contigs(contig_file,filter_params,type_reads,num_threads,*reads):
    filtered_contig_file=contig_file[:contig_file.rfind('/')+1]+'filtered_'+contig_file[contig_file.rfind('/')+1:]
    shutil.copy(contig_file,filtered_contig_file) #makes a copy of the file renamed with "filtered_"
    print('\nFiltering reads...')
    params=filter_params.split('|')
    for i in params:
        code,val=i.split('=')
        if code=='m':
            print('... Filter out contigs < '+val+' nucleotides in length.')
            filtered_contig_file=filter_contigs_by_size(filtered_contig_file,'>',val)
            print('\t'+str(get_number_of_sequences(filtered_contig_file))+' remain.')
        elif code=='M':
            print('... Filter out contigs > '+val+' nucleotides in length.')
            filtered_contig_file=filter_contigs_by_size(filtered_contig_file,'<',val)
            print('\t'+str(get_number_of_sequences(filtered_contig_file))+' remain.')
        elif code=='c':
            print('... Filter out contigs with coverage < '+val+'x.')
            filtered_contig_file=filter_by_coverage(filtered_contig_file,val,type_reads,reads[0][0],reads[0][1])
            print('\t'+str(get_number_of_sequences(filtered_contig_file))+' remain.')
        elif code=='cov':
            print('... Filter out contigs with a SPAdes cov value < '+val+'x.')
            filtered_contig_file=filter_by_spades_cov(filtered_contig_file,val)
            print('\t'+str(get_number_of_sequences(filtered_contig_file))+' remain.')
        else:
            print('... Filter out contigs that do not show homology to the sequences of interest. File: '+val+'.')
            filtered_contig_file=filter_by_seq_of_interest(filtered_contig_file,val,num_threads)
            print('\t'+str(get_number_of_sequences(filtered_contig_file))+' remain.')
    return filtered_contig_file

#CHECKS THAT THE OUTPUT PATH EXISTS. IF IT ALREADY EXISTS, IT MAKES ANOTHER ONE.
def check_output_path(output_path):
    #check that output_path exists; if not, make it
    if not os.path.exists(output_path):         ## does not exist
        os.makedirs(output_path+'/temp')
        if not os.path.exists(output_path+'/temp'):
            print("Cannot write to output folder. Please check the output path and permissions.")
            return False,''
    return True,output_path

#rename contigs to consistent numbering after processing
def rename_contigs(contig_file,output_path):
    contigs=[]
    contigs=list(SeqIO.parse(contig_file,'fasta'))
    outputf=open(output_path+'/final_contigs.fasta','w')
    counter=1
    for i in contigs:
        i.name=''
        i.description=''
        i.id='contig_'+str(counter)
        counter=counter+1
    SeqIO.write(contigs,output_path+'/final_contigs.fasta','fasta')
    del contigs[:]
    return output_path+'/final_contigs.fasta'

#RUN PIPELINE
def run_virMine(assembly,contig_file,type_reads,reads,filter_params,output_path,bact_db,viral_db,num_threads,fast):
    if assembly!='provided':
        if qc_and_assembly(assembly,str(num_threads),type_reads,reads)==False:
            print('\nError in QC and assembly. virMine did not run.')
            log('Error in QC and assembly. virMine did not run.',output_path+'/virMine_log.txt')
            return False
        else:
            if assembly=='megahit':
                contig_file=output_path+'/temp/assembly/final.contigs.fa'
            elif assembly=='spades':
                contig_file=output_path+'/temp/assembly/scaffolds.fasta'
            elif assembly=='metaspades':
                contig_file=output_path+'/temp/assembly/scaffolds.fasta'
            elif assembly=='all3':
                contig_file=output_path+'/temp/assembly/scaffolds.fasta'
    else: #assembly is provided, copy it to the output_path
        temp_output_file=output_path+'/temp/'+contig_file[contig_file.rfind("/")+1:]
        shutil.copyfile(contig_file,temp_output_file)
        contig_file=temp_output_file
        log('Assmebly provided & copied to: '+contig_file,output_path+'/virMine_log.txt')

    num_contigs_from_assembly=get_number_of_sequences(contig_file)
    log('Number of contigs in assembly file '+contig_file+': '+str(num_contigs_from_assembly)+'.',output_path+'/virMine_log.txt')
    if num_contigs_from_assembly==0:
        print('Error with assembly.')
        return False
    contig_file=qc_contig_file(contig_file) #reformat contigs file
    if contig_file==False:
        print('\nError with QC of contig file. Analysis cannot continue.')
        return False
    print('Assembly produced '+str(num_contigs_from_assembly)+' contigs')


    #RUN FILTERS
    if filter_params!='':
        log('Running Filters...',output_path+'/virMine_log.txt')##
        contig_file=filter_contigs(contig_file,filter_params,type_reads,num_threads,reads)
        num_contigs_after_filtering=get_number_of_sequences(contig_file)
        print('Filtering resulted in '+str(num_contigs_after_filtering)+' contigs')
        log('Filtering resulted in '+str(num_contigs_after_filtering)+' contigs.',output_path+'/virMine_log.txt')
        if num_contigs_from_assembly==num_contigs_after_filtering:
            print('Filters did not remove any contigs.')
            log('Filters did not remove any contigs.',output_path+'/virMine_log.txt')
        else:
            print('Filters removed '+str(num_contigs_from_assembly-num_contigs_after_filtering)+' contigs.')
            log('Filters removed '+str(num_contigs_from_assembly-num_contigs_after_filtering)+' contigs.',output_path+'/virMine_log.txt')
    else:
        print('No filters applied.')
        log('No filters applied.',output_path+'/virMine_log.txt')##

    #copy filtered reads to the main folder
    contig_file=rename_contigs(contig_file,output_path)
        
    ##RUN GLIMMER AND PREDICTIONS
    orf_file=run_glimmer(contig_file,output_path)
    if orf_file==False:#there was an error in predicting genes/ no genes were predicted
        num_predicted_orfs=0
    else:
        num_predicted_orfs=get_number_of_sequences(orf_file)
    if num_predicted_orfs==0:
        print('No ORFs were predicted.')
        log('No ORFs were predicted.',output_path+'/virMine_log.txt')##
        return False
    else:
        log('# Predicted ORFs: '+str(num_predicted_orfs),output_path+'/virMine_log.txt')

    if fast=='False':
        num_called_bacteria, num_called_viral, num_called_unkn, v_contigs, unkn_contigs = determine_contig_scores(contig_file,orf_file,bact_db,viral_db,output_path,num_threads)
        write_out_blast_results(output_path+'/all_viral.blastx',output_path+'/viral_contigs.blastx',v_contigs)
    else:
        num_called_bacteria, num_called_viral, num_called_unkn, v_contigs, unkn_contigs = determine_contig_scores_p(contig_file,orf_file,bact_db,viral_db,output_path,num_threads)
        write_out_blast_results(output_path+'/all_viral.blastp',output_path+'/viral_contigs.blastp',v_contigs)
    log('\nNumber of Contigs Classified:\nNumber Bacterial Contigs: '+str(num_called_bacteria)+'\nNumber Viral Contigs: '+str(num_called_viral)+'\nNumber Unknown Contigs: '+str(num_called_unkn),output_path+'/virMine_log.txt')

    write_out_contigs_of_interest(contig_file,output_path+'/viral_contigs.fasta',v_contigs)
    log('Viral contigs written to: '+output_path+'/viral_contigs.fasta',output_path+'/virMine_log.txt')
    write_out_contigs_of_interest(contig_file,output_path+'/unkn_contigs.fasta',unkn_contigs)
    log('Unknown contigs written to: '+output_path+'/unkn_contigs.fasta',output_path+'/virMine_log.txt')
    write_out_orfs_for_contigs_of_interest(orf_file,output_path+'/viral_contig_orfs_nt.fasta',v_contigs)
    log('ORFs (nt) for viral contigs written to: '+output_path+'/viral_contig_orfs_nt.fasta',output_path+'/virMine_log.txt')
    final_orfs=list(SeqIO.parse(output_path+'/viral_contig_orfs_nt.fasta','fasta'))
    orf_output=open(output_path+'/viral_contig_orfs_aa.fasta','w')
    for i in final_orfs:
        orf_output.write('>'+str(i.id)+'\n')
        orf_output.write(str(i.seq.translate())+'\n')
    orf_output.close()
    log('ORFs (aa) for viral contigs written to: '+output_path+'/viral_contig_orfs_nt.fasta',output_path+'/virMine_log.txt')##
    return True


ts=time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
command_message='Running virMine.\nRun started at '+st+'\n\nCommand:\n'
command_message+='python'
for i in sys.argv:
    command_message+=' '+i


#SET VARIABLES
if args.assembler is None:
    assembly='provided'
    contig_file=args.assembled_contigs
    type_reads=None
    reads=None
    if args.single_reads is not None:
        type_reads='se'
        reads=args.single_reads
    else:
        type_reads='pe'
        reads=args.paired_end_reads
else:
    contig_file=''
    assembly=args.assembler
    if args.single_reads is not None:
        type_reads='se'
        reads=args.single_reads
    else:
        type_reads='pe'
        reads=args.paired_end_reads
output_path=args.output_path

if args.blast_option is None:
    fast=True
else:
    if args.blast_option=='blastx':
        fast=False
    else:
        fast=True

if args.num_threads is None:
    num_threads=1
else:
    num_threads=args.num_threads
if num_threads>multiprocessing.cpu_count():
    num_threads=1

filter_params=''
if args.min_contig_size is not None:
    filter_params+='m='+str(args.min_contig_size)
if args.max_contig_size is not None:
    if len(filter_params)>0:
        filter_params+='|'
    filter_params+='M='+str(args.max_contig_size)
if args.min_coverage is not None:
    if assembly=='provided':
        if type_reads is None:
            print('Cannot filter by coverage. No reads were provided.\n')
        else:
            if len(filter_params)>0:
                filter_params+='|'
            filter_params+='c='+str(args.min_coverage)
    else:
        if len(filter_params)>0:
            filter_params+='|'
        filter_params+='c='+str(args.min_coverage)        
if args.min_SPAdes_cov is not None:
    if assembly=='provided':
        print('Filtering on SPAdes cov values not available for assemblies provided by the user.\n')
    else:
        if len(filter_params)>0:
            filter_params+='|'
        filter_params+='cov='+str(args.min_SPAdes_cov)
if args.genes_of_interest is not None:
    if len(filter_params)>0:
        filter_params+='|'
    filter_params+='g='+str(args.genes_of_interest)
bact_db_seqs=args.non_viral
viral_db_seqs=args.viral


#CHECK PATH, FILES, and PATHS PROVIDED BY THE USER
x=check_installs() #check that the PATH variables are set correctly
if x==False:
    print('Path variables must be set prior to execution.\nCheck PATH and try again.')
else:
    #check reads file (if starting from raw data); check assembly if starting from there
    if assembly=='provided':
        fileQC,errmsg=check_files(contig_file)
    else:
        if type_reads=='se':
            fileQC,errmsg=check_files(reads)
        else:
            fileQC,errmsg=check_files(reads[0],reads[1])
    if fileQC==False:
        print(errmsg)
    else:
        #check DB locations
        fileQC,errmsg=check_files(bact_db_seqs,viral_db_seqs)
        if fileQC==False:
            print('DATABASE SEQUENCE ERROR: '+errmsg)
        else:
            x,output_path=check_output_path(output_path)
            log('Output will be written to: '+output_path+'.',output_path+'/virMine_log.txt')
            if x==True:
                #make DBs
                bact_db=output_path+'/temp/'+'nonviral'
                makeblast_command='makeblastdb -in '+bact_db_seqs+' -out '+bact_db+' -title '+bact_db+' -dbtype prot'
                os.system(makeblast_command)          
                viral_db=output_path+'/temp/'+'viral'
                makeblast_command='makeblastdb -in '+viral_db_seqs+' -out '+viral_db+' -title '+viral_db+' -dbtype prot'
                os.system(makeblast_command)
                fileQC,errmsg=check_files(bact_db+'.pin',viral_db+'.pin')
                if fileQC==False:
                    print('DATABASE ERROR: '+errmsg)
                    log('DATABASE ERROR: '+errmsg,output_path+'/virMine_log.txt')
                else:
                    log(command_message+'\n\n*****************************************\n\n',output_path+'/virMine_log.txt')
                    log('Viral & nonviral BLAST databases created:\n\t'+bact_db+'\n\t'+viral_db,output_path+'/virMine_log.txt')
                    runOK=run_virMine(assembly,contig_file,type_reads,reads,filter_params,output_path,bact_db,viral_db,num_threads,fast)
                    if runOK==True:
                        print('Your job has finished running.\n')
                        log('virMine completed.',output_path+'/virMine_log.txt')
                    else:
                        print('Error occurred. Reference virMine_log.txt file.')
                        log('Error occurred. Reference virMine_log.txt file.',output_path+'/virMine_log.txt')
                    ts=time.time()
                    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
                    log('virMine completed at: '+st,output_path+'/virMine_log.txt')
