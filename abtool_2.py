#!/usr/local/bin/python

########################################################################
# Antibody pairing and fasta file generator
#
# Written by: Edwin Kamau 
# Last modified: 3.15.2019 by Edwin Kamau
#
########################################################################
# import all the necessary functions required to make functional calls

import sys
import os
import re
import subprocess
import platform
import glob 

if len(sys.argv)<7:
    print("Usage: abtool.py  <path to json_file>  <outname> <organism genus or species name> <organism common_name> <receptor_type> <threshold> <version>\n\n json_file: NGS sequencing output file in JSON format eg from 10x Genomics platform\n Outname:prefix that will be used to name the output files. Defaults to output\n Organism genus or species name: generic or specific name of the organism of interest eg. Homo or sapiens for human Mus or musculus for mouse. See IMGT naming convention for guidance\n Common_name:organism common name e.g human,mouse,...\n Very important as igblast modules require proper naming convention\n Receptor:Type of receptor query either immunoglobulin(IG) or T cell receptor (TR)\n Threshold: will output paired clones whose heavy chain percent identity is below the indicated value\n Version:Optional input the for prefered igblast to use version.Program installs version 1.8.0\n")
    sys.exit(0)

# GRAB SYSTEM LEVEL VARIABLES
path=sys.argv[1]
outname=sys.argv[2]
working_dir=os.getcwd()
set1='imgt'
sp=sys.argv[3] # Homo sapiens
cn=sys.argv[4] # common name human
receptor=sys.argv[5] # 'IG/TR'
threshold = sys.argv[6]
release='1.8.0'

def output(dataframe):
    '''
    This will be output of the created dataframe 
    '''
    # load modules 
    from collections import Counter
    import pandas as pd

    # make a copy of the dataframe
    df=dataframe.copy()
    df=df.set_index('clonotype') # change index to 'clonotype'
    df.to_csv(outname + '.csv')
    print('Generating final output based on threshold given')

    # final columns to output f
    reduced=['contig_name_hc','contig_name_lc','v_gene_hc','v_gene_lc','v_gene_percent_identity_hc','v_gene_percent_identity_lc','cdrs_hc','cdrs_lc','d_gene_hc', 'j_gene_hc','j_gene_lc', 'isotype_hc']

    all=['contig_name_hc','contig_name_lc','v_gene_hc','v_gene_lc','v_gene_percent_identity_hc','v_gene_percent_identity_lc','cdrs_hc',
        'cdrs_lc','d_gene_hc', 'j_gene_hc','j_gene_lc', 'isotype_hc',
        'cdr3_imgt_aa_hc','cdr3_imgt_aa_lc','sequence_hc','sequence_lc']

    # Identify paired clones and pair them 
    pairs=[key for key,value in Counter(df.index).items() if value>1]
    paired=df[df.index.isin(pairs) & df.chain.isin(['VH'])].join(df[df.index.isin(pairs) & df.chain.isin(['VL','VK'])],how='outer', rsuffix='_lc',lsuffix='_hc').dropna(subset=['contig_name_hc','contig_name_lc'])[all]
    paired.to_csv(outname + 'all_paired.csv')
    out=paired[(paired.v_gene_percent_identity_hc<threshold)][fields_output].sort_values(by=reduced[4])

    out.to_csv(outname + 'by_percent.csv')
    #multiple_pair=[key for key,value in Counter(df.index).items() if value>2]
    #multiple_paired=df[df.index.isin(multiple_pair) & df.chain.isin(['IGH'])].join(df[df.index.isin(multiple_pair) & df.chain.isin(['IGK','IGL'])],how='outer',rsuffix='_lc').join(df[df.index.isin(multiple_pair) & df.chain.isin(['IGK','IGL'])],how='outer',rsuffix='_lc2')
            # remove duplicated rows
     #       multiple_paired=multiple_paired[multiple_paired.IGHV_lc !=multiple_paired.IGHV_lc2].drop_duplicates('cdr3')

def create_fasta_file(path,outname='output'):
    '''
    function requires file to be analyzed. Assumes a particular annotation format thus may need to be modified if the format changes 
    if outname is not given it will use output as the prefix
    outputs: DNA and protein fasta sequences as single or multiple paired
             single and multiple paired chains as csv format  
    '''
    try:
    	import pandas as pd
    except Exception as e:
      print(e)
      setup()
    try:
        if path:
            print('Begin processing %s'% path)
            # load file
            db=pd.read_json(path,orient='records')
            # REMOVE ALL NON-PAIRED chains
            df=db.copy()
            df=df[df.duplicated('clonotype',keep=False)]# remove non paired clones
            # get chain type, gene name and gene start and end from annotation
            f=open(outname +'_dna_seq.fasta','w')
            for rec in df.index:
                sequence=str(df.loc[rec]['sequence']) 
                df.loc[rec,'chain']=df.loc[rec,'annotations'][0]['feature']['chain']# 
                gene_id=[d['feature']['gene_name'] for d in df.loc[rec,'annotations']]# grab gene names 
                position=[(s['contig_match_start'],s['contig_match_end']) for s in df.loc[rec,'annotations']]
                combine={} # Combine gene names and start and end and filter 5' UTR 
                for ids in range(len(gene_id)):
                    if gene_id[ids] not in combine:
                        combine[gene_id[ids]]=position[ids]
                for key,value in combine.items():
                    if 'V' in key:
#                        df.loc[rec,'IGHV']=key
                        df.loc[rec,'v_start']=value[0]
                        df.loc[rec,'sequence']=sequence[value[0]:]# remove 5' UTRdf=df.set_index('clonotype') # change index to 'clonotype'
                        df.loc[rec,'v_end']=value[1]
                    #if 'D' in key and len(key)>4:
 #                       df.loc[rec,'IGHD']=key
                        #df.loc[rec,'d_start']=value[0]
                        #df.loc[rec,'d_end']=value[1]
                    if 'J' in key:
 #                       df.loc[rec,'IGHJ']=key
                        df.loc[rec,'j_start']=value[0]
                        df.loc[rec,'j_end']=value[1]
                    if len(key)==5 and key[3]!='J' or len(key)==4:
                        df.loc[rec,'isotype']=key
                        df.loc[rec,'iso_start']=value[0]
                        df.loc[rec,'iso_end']=value[1]
                        f.write('>%s|%s|%s|%s\n%s\n' %(df.loc[rec,'contig_name'],
                              df.loc[rec,'isotype'],df.loc[rec,'v_start'],
                              df.loc[rec,'j_end'],df.loc[rec,'sequence']))
            f.close()
            print("Done processing %s\n" % (path))
    except Exception as e:
        print(e)
        sys.exit()

def igblast(seqfile,igblastfile):
    '''
    Summarise each Igblast analysis on a single tab-separated line.',
    epilog='The program creates two output files: <tag>_n.fasta (nucleotide analysis) and <tag>_n.fasta (amino acid analysis)\nRequires fasta file output of abtool
    fasta file and runigblast produced with outfmt 3 and a prefix for output name )
    parser.add_argument('seqfile', help='file containing the sequences analysed by IgBlast (FASTA)')
    '''
    # LOAD REQUIRED MODULES

    from Bio import SeqIO
    from Bio.Seq import Seq
    import traceback
    import pandas as pd 

    fieldorder = ['fr1', 'cdr1', 'fr2', 'cdr2', 'fr3', 'cdr3']
    fieldnames = ['clonotype','contig_name','chain','v_gene', 'v_gene_percent_identity','d_gene', 'j_gene', 'isotype','cdrs','cdr3_imgt', 'Functionality', 'Stop_Codon', 'V_J_frame','Strand', 'fr1_imgt', 'cdr1_imgt', 'fr2_imgt', 'cdr2_imgt', 'fr3_imgt', 'junction_imgt', 'Notes','fr1_imgt_aa', 'cdr1_imgt_aa', 'fr2_imgt_aa', 'cdr2_imgt_aa', 'fr3_imgt_aa', 'cdr3_imgt_aa','sequence']

    seqs = {}
    df=[]
    print('Performing igblast analysis on %s output' %(igblastfile))
    for seq_record in SeqIO.parse(seqfile, 'fasta'):
        clone=seq_record.id.split('_')[0]
        seqs[clone] = str(seq_record.seq.upper())

    with open(igblastfile, 'r') as fi:
        res = {}
        junc = {}
        lastline = ''
        alignments = False
        
        for line in fi:
            try:
                rawline = line
                line = line.strip()
                items = line.split('\t')
                
                if 'Query= ' in line:                        
                    contig,isotype,v_start,j_end=(line.replace('Query= ', '').strip()).split('|')[0:]
                    res['clonotype'] =contig.split('_')[0]
                    res['contig_name']=contig
                    res['isotype']=isotype
                    res['v_start']=v_start
                    res['j_end']=j_end
                    res['Notes'] = ''
                    alignments = False
                elif 'V-(D)-J rearrangement summary for query sequence' in lastline:               
                    res['v_gene'] = items[0]

                    if 'D gene' in lastline:
                        res['d_gene'] = items[1]
                        bump = 1
                    else:
                        res['d_gene'] = ''
                        bump = 0
                        
                    res['j_gene'] = items[bump+1]
                    res['chain'] = items[bump+2]
                    res['Stop_Codon'] = items[bump+3]
                    res['V_J_frame'] = items[bump+4]                
                    res['Functionality'] = 'productive' if items[bump+5] == 'Yes' else 'unproductive'
                    res['Strand'] = items[bump+6]
    
                    if ',' in res['v_gene']:
                        selected_vgene = res['v_gene'].split(',')[0]
                    else:
                        selected_vgene = res['v_gene']
    
                    if ',' in res['d_gene']:
                        selected_dgene = res['d_gene'].split(',')[0]
                    else:
                        selected_dgene = res['d_gene']
    
                    if ',' in res['j_gene']:
                        selected_jgene = res['j_gene'].split(',')[0]
                    else:
                        selected_jgene = res['j_gene']
    
                elif 'V-(D)-J junction details' in lastline:               
                    junc['clonotype'] = res['clonotype']
                    junc['Functionality'] = res['Functionality']
                    
                    for i in range(0, len(items)):
                        if items[i] == 'N/A':
                            items[i] = ''
                    
                    if 'D region' in lastline:
                        junc['3\'v_region'] = items[0]
                        junc['n_region'] = (items[1].replace('(', '')).replace(')', '')
                        junc['d_region'] = items[2]
                        junc['n2_region'] = (items[3].replace('(', '')).replace(')', '')
                        junc['5\'j_region'] = items[4]
                    else:
                        junc['3\'v_region'] = items[0]
                        junc['n_region'] = (items[1].replace('(', '')).replace(')', '')
                        junc['d_region'] = ''
                        junc['n2_region'] = ''
                        junc['5\'j_region'] = items[2]

                elif 'Total' in line and 'N/A' in line:
                    res['v_gene_percent_identity']=items[-1]                       
                elif 'Alignments' in line:
                    alignments = True
                    al_narrative = ''
                    al_query_seq = ''
                    al_query_start = 0
                    al_query_frag = ''
                    al_query_frag_start = 0
                    al_query_j_start = 0
                    al_query_end = 0
                    al_vmatch_start = 0
                    al_jmatch_start = 0
                    al_jmatch_end = 0
                    raw_fr3 = ''
                    
                    if 'j_gene' in res and res['j_gene'] != 'N/A':
                        jmatch = res['j_gene'] if len(res['j_gene']) == 1 else res['j_gene'].split(',')[0]
                    else:
                        jmatch = ''
                        
                elif alignments:
                    items = rawline.split()
                    if len(items) > 1:
                        if 'Query_' in line and len(items) == 4:
                            if len(al_query_seq) == 0:
                                al_query_start = int(items[1])
                            al_query_frag_start = int(items[1])
                            al_query_frag = items[2]
                            al_query_seq += items[2]
                            al_query_end = int(items[3])
                            
                            p = rawline.find(items[2])
                            al_narrative += lastline[p:p + len(items[2])]
                            
                        elif len(items) == 7 and jmatch == items[3]:
                            al_jmatch_end = int(items[6])
                            if al_jmatch_start == 0:
                                al_jmatch_start = int(items[4])
                            if al_query_j_start == 0:
                                query_frag_gaps = 0
                                for i in range(0, len(items[5])):
                                    if al_query_frag[i] == '-':
                                        query_frag_gaps += 1
                                    if items[5][i] != '-':
                                        al_query_j_start = al_query_frag_start + i - query_frag_gaps
                                        break                    
                        elif 'Lambda' in line:
                            alignments = False
                            if len(al_query_seq) != len(al_narrative) and len(al_narrative) > 0 and len(al_query_seq) > 0:
                                print('Error: alignment sequence and narrative are out of synch in id %s' % res['clonotype'])
                            
                            # field labels may be truncated if the fields are short. In the worst case, this can be right down
                            # to a single letter: <C> or <F>. I have not seen any instances of anything shorter than that, but
                            # we should check.
                            
                            def findextents(s_needle, e_needle, haystack):
                                extents = []
                                i = 0
                                haystack_len = len(haystack)
                                while i < haystack_len:
                                    i = haystack.find(s_needle, i)
                                    if i < 0:
                                        return extents
                                    j = haystack.find(e_needle, i + len(s_needle)) + 1
                                    if j < 0:
                                
                                        print("Broken field delimiter in %s: no closing %s." % (haystack, e_needle))
                                        return extents
                                    extents.append((i, j))
                                    i = j
                            
                            extents = findextents('<', '>', al_narrative)
                            
                            if extents is not None:
                                field_values = []
                                field_indeces = []
                                for (start, end) in extents:
                                    fdesc = al_narrative[start:end]
                                    field_values.append(al_query_seq[start:end])
                                    if '<' in fdesc[1:]:
                                        print ("Broken field delimiter in %s: unexpected '<" % al_narrative)
                                        break
                                    for i, v in enumerate(fieldorder):
                                        fval = -1
                                        if v.upper() in fdesc:
                                            fval = i
                                            break
                                    field_indeces.append(fval)
                                # we require at least one unambiguous field
                                start_ind = -1
                                for i, v in enumerate(field_indeces):
                                    if v >= 0:
                                        start_ind = v - i
                                        break
                                        
                                fields_correct = True
                                if start_ind < 0:
                                    print ("%s: Alignment parser can't infer any fields in '%s'" % (res['clonotype'], al_narrative))
                                else:
                                    for i,v in enumerate(field_indeces):
                                        if v < 0:
                                            field_indeces[i] = start_ind + i
                                        elif v != start_ind + i:
                                            print("%s: Alignment fields out of order, or missing field, in %s" % (res['clonotype'], al_narrative))
                                            fields_correct = False
                                            break
                                # group sequence by framework or cdr and translate 
                                # will ignore clones with missing framework or cdr regions
                                if fields_correct and len(field_indeces)==6: 
                                    for ind in field_indeces:
                                        res[fieldorder[ind] + '_imgt'] = field_values[ind].replace('-', '')
                                        res[fieldorder[ind] + '_imgt_aa'] = str(Seq(field_values[ind].replace('-', '')).translate())
                                        if fieldorder[ind] == 'cdr3':
                                            (start, end) = extents[ind]
                                            if len(al_query_seq) > end+3:
                                                res['junction_imgt'] = al_query_seq[start-3:end+3].replace('-', '')
                                            
                                
                                    if res['clonotype'] not in seqs:
                                        raise ValueError("Sequence %s not found in FASTA file" % res['clonotype']) 
                                        
                                    seq = seqs[res['clonotype']]
                                    res['sequence'] = seq
                                    res['cdrs']=str(len(res['cdr1_imgt_aa']))+ '.' + str(len(res['cdr2_imgt_aa']))+'.'+str(len(res['cdr3_imgt_aa']))
                                    df.append(res.copy())
                                
                            # Sometimes IgBlast fails to match the entire fr1 region.
                            # Extend fr1 back to cover as much as possible of the v-gene.
                            # Remember indeces reported by IgBlast are 1-based not 0-based.
                            if 'fr1_imgt' in res and al_query_start > 1 and al_vmatch_start > 1:
                                extent = min(al_query_start, al_vmatch_start)
                                res['fr1_imgt'] = seq[al_query_start - extent - 1: al_query_start - 1] + res['fr1_imgt']
                                while (len(res['fr1_imgt']) % 3) != 0:
                                    res['fr1_imgt'] = res['fr1_imgt'][1:]
                                    
                            if 'junction_imgt' in res:
                                inner_junc = junc['n_region'] + junc['d_region'] + junc['n2_region']
                                if inner_junc in res['junction_imgt']:
                                    p = res['junction_imgt'].find(inner_junc)
                                    junc['3\'v_region'] = res['junction_imgt'][:p] if p > 0 else ''
                                    junc['5\'j_region'] = res['junction_imgt'][p + len(inner_junc):] if p + len(inner_junc) < len(res['junction_imgt']) else ''
                                    
                                    if junc['3\'v_region'] + inner_junc + junc['5\'j_region'] != res['junction_imgt']:
                                        print("%s: junction misalignment" % res['clonotype'] ) 
                                        res['Notes'] += 'Error: junction misalignment'
                                else:
                                    # If IgBLAST does not find a V-gene alignment which extends as far as the first Cys of the junction, it creates an N-region which covers
                                    # the first Cysteine and extends downstream beyond the junction. As the junction analysis is inaccurate, we mark such (rare) occurrences 
                                    # as non-productive, even though IgBLAST was able to determine the junction in the alignment.
                                    print("%s: junction analysis %s does not match inferred junction %s" % (res['clonotype'], inner_junc, res['junction_imgt']))
                                    res['Notes'] += ' Junction analysis does not match inferred junction. Possibly first Cysteine was not identified in V-gene alignment.'
                                    res['Functionality'] = 'unproductive'
                            else:
                                res['Functionality'] = 'unproductive'                        
                                if jmatch == '':
                                    res['Notes'] += 'J-gene not identified. '
                                if 'fr3_imgt' not in res:
                                    res['Notes'] += 'fr3 not identified. '
                        
                            if len(res) > 0:
                                for k,v in res.items():
                                    if '_length' in k:continue
                                    if 'N/A' in v:
                                        res[k] = ''         # blank out any N/As from IgBLAST for compatibility with imgt
        
                                junc['clonotype'] = res['clonotype']
                                junc['Functionality'] = res['Functionality']
                                junc['v_gene'] = res['v_gene']
                                junc['d_gene'] = res['d_gene']
                                junc['j_gene'] = res['j_gene']
                                junc['Notes'] = res['Notes']
                                
                                if 'junction_imgt' in res:
                                    junc['junction'] = res['junction_imgt']
                                
                                for k,v in junc.items():
                                    if 'N/A' in v:
                                        junc[k] = ''         # blank out any N/As from IgBLAST for compatibility with imgt
                                
                lastline = rawline
            
            except:
                id = res['clonotype'] if 'clonotype' in res else '<unknown>'
                print("Error parsing sequence %s:\n%s" % (id, traceback.format_exc()))
    df=pd.DataFrame.from_records(df)

    output(df[fieldnames])

def runigblast(igblastcmd,germpath,optional):
    '''
    generate igblast file 
    '''
    file_root='imgt'+ '_' + cn + '_' + receptor
    vfile = germpath + '/' + file_root + 'V'
    dfile = germpath + '/' + file_root + 'D'
    jfile = germpath + '/' + file_root + 'J'
    seqfile = glob.glob(working_dir + '/' + outname + '*dna*.fasta')[0]

    print('Begin mapping %s to imgt %s database' % (seqfile,cn))

    logfile=outname + '_igblast.log'

    cmd = '%s -germline_db_V %s -germline_db_D %s -germline_db_J %s -organism %s -domain_system %s -query %s -auxiliary_data %s/%s_gl.aux -outfmt 3 >%s' % (igblastcmd, vfile, dfile, jfile, cn, set1, seqfile, optional,cn, logfile)
    print('\nRunning igblast\n%s' % cmd)
    subprocess.call(cmd, shell=True)

    print('Finished mapping sequences to imgt')
    igblast(seqfile,logfile)

def imgt(igblastpath,igblastcmd,folder):
    '''
    downloads imgt germline sequences from imgt and makes database base on specified organism
    '''
    # make germline folder path
    try:
        germpath=igblastpath + '/%s' % set1
        def make_db():
            recs = list(SeqIO.parse('%s/imgt_germlines.fasta' % germpath, 'fasta'))
            fastafiles = []
            outrecs = {}
            for rec in recs:
                d = rec.description.split('|')
                if len(d) > 4 and receptor in d[1] and sp in d[2] and d[4] in ['V-REGION', 'D-REGION', 'J-REGION']:
                    seg = d[4][0]
                    rec.description = d[1]
                    rec.id = d[1]
                    key = d[1][:2] + seg
                    if key not in outrecs:
                        outrecs[key] = []
                    outrecs[key].append(rec)
            for fn, r in outrecs.items():
                fastafile = '%s/imgt_'% germpath + cn + '_' + fn + '.fasta'
                SeqIO.write(r, fastafile, 'fasta')
                fastafiles.append(fastafile)
            for fn in fastafiles:
                print('Processing germline file %s' % fn)
                cmd = '%s/makeblastdb -parse_seqids -dbtype nucl -in %s -out %s' % (igblastpath, fn, fn.replace('.fasta', ''))
                subprocess.call(cmd, shell=True)   
            print('Germline file processing complete.')
            #subprocess.call('touch %s/complete' % path, shell=True)

        if os.path.isfile(germpath + '/%s_germlines.fasta' %(set1)):
            if os.path.isfile(germpath + '/%s_%s_%sV.nhr' %(set1,cn,receptor)):
                print('imgt_germline file previously downloaded and processed for species : %s' % cn)
            if not os.path.isfile(germpath + '/%s_%s_%sV.nhr' %(set1,cn,receptor)):
                print('making %s %s database for' % (cn,set1))
                make_db()

        if not os.path.exists(germpath):
            os.makedirs(germpath)
            print ('Installing germline files in  %s.' % germpath)

            subprocess.call("wget -P %s -O %s/imgt_germlines.fasta http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP" % (germpath,germpath), shell=True)
            make_db()
    except exception as e:
        print(e)
    abtool(path,outname)
    runigblast(igblastcmd,germpath,folder)

def setup():
    '''
    sets up igblast executables in home director and the required files internal_data, optional files and imgt database in $HOME/igblast
    '''
    # dependent files generator
    print("Setting up your system for igblast executables and dependencies\nExecutables installed in home directory")

    def dependent(igblastpath,folder):
        '''
        creates dependent files for igblast 
        '''
        path=igblastpath +'/' + folder
        if not os.path.exists(path):# check if directory exist
            os.makedirs(path) # create directory if non existent
            print ('Igblast executables installed but no %s dependency files\nInstalling the dependencies' % folder)
            subprocess.call("wget -P %s -nH --cut-dirs=5 -r ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/%s/ " %(path,folder),shell=True)#download dependency files
        if 'internal_data' in path:# create a link to internal data doesnt work if this is not done
            print('linking %s to %s as %s' % (path,working_dir,folder))
            subprocess.call('ln -s %s %s/%s' % (path,working_dir,folder),shell=True) # download dependency file
            print ('Finished installing %s files' % folder)
        if 'optional_file' in path:# return optional file 
            print(path)
            return path

    home=os.environ['HOME']


    igblastbin=home + '/ncbi-igblast-' + release + '/bin/' # path of the executables

    igblastcmd=igblastpath + '/igblastn' # executable file

    if os.path.exists(igblastpath):
        if os.path.isfile(igblastcmd): # is executable installed?
            if os.path.exists(igblastpath + '/internal_data'): # is dependent file installed
                if os.path.exists(igblastpath + '/optional_file'): # is optional file installed
                    print('igblast executables and dependent files previously installed')
                if not os.path.exists(igblastpath + '/optional_file'):
                    folder=dependent(igblastpath,'optional_file')
            if not os.path.exists(igblastpath + '/internal_data'): # are dependent file installed
                dependent(igblastpath,'internal_data')

    if not os.path.exists(igblastpath): # no igblast executables
        print('Installing igblast executables and dependent files in the %s directory' % home)

        url= 'ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/%s/ncbi-igblast-%s-x64-%s.tar.gz' %(release,release,system)

        extension=url.split('/')[-1]   

        subprocess.call("wget -P %s %s ; tar -xf %s/%s -C %s" %(home,url,home,extension,home),shell=True)
        for depend in ['internal_data','optional_file']:
            folder=dependent(igblastpath,depend)

    imgt(igblastpath,igblastcmd,folder) # pass to imgt to make db a

def install():
    '''
    This function will install the required packages required for processing input file and perform other needed computational processes.
    write now its works for macos and linux system. 
    '''
    try:
        if platform.system()=='Darwin':
            system='macosx'
        if platform.system()=='Linux':
            system='linux'
        if not system:
            print('Your platform %s is not supported for installation of packages and modules\nPlease check if you have wget,python3, and python packages\n scipy,numpy,pandas,and biopython and either install or rerun the program %s\n' %(platform.system(), platform.system())
            sys.exit()
        # define brew command
        brew = '$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)'
        if system =='macosx': 
        # check and install dependencies
            for prog in ['xcode-select','brew','python3','wget']:
                if subprocess.call('which %s' % prog,shell=True)==0:
                    print('%s already installed'% prog)
                    if prog == 'brew':
                        print('Updating %s please wait\n' % prog)
                        subprocess.call('%s update' % prog,shell=True)

                if subprocess.call('which %s' % prog,shell=True)!=0:
                    try:    
                        print('%s is not installed in your system')
                        if prog == 'xcode-select':
                                print('Installing %s\will require permission to install\nplease respond to dialog boxes that will popup\n Password: The password used to login into the machine otherwise you should get the machine owner/institution and rerun the script to install %s' %(prog,prog))
                                subprocess.call('%s --install',shell=True)
                                print('%s installed' % prog )
                        if prog =='brew':
                                print('installing %s' % prog)
                                subprocess.call('ruby -e %s' % brew, shell=True)
                                print('brew makes app installation easier\nUse it to install new packages or upgrade downloadedones\n')
                        if prog.startswith('python3'):
                            print('installing %s' % prog)
                            if subprocess.call('which brew',shell=true)==0:
                                subprocess.call('brew install %s' % prog, shell=True) # should install pip3 as well
                            # install python dependent packages
                            if subprocess.call('which pip3',shell=True)==0:
                                subprocess.call('pip3 install scipy numpy pandas biopython',shell=True)
                        if prog.startswith('wget'):
                            print('installing %s' % prog)
                            if subprocess.call('which brew',shell=true)==0:
                                subprocess.call('brew install %s' % prog, shell=True)
                    except Exception as e:
                        print(e)
                        print('Some or all the packages did not install correctly and require manual installation\nFollow the instructions of the popup web page to install them\nThen run the program again')
                        subprocess.call('open https://www.howtogeek.com/211541/homebrew-for-os-x-easily-installs-desktop-apps-and-terminal-utilities/')
                        sys.exit()

    except Exception as e:
        print(e)
    print('Your system has the required modules to run this program\nPlease update the modules regularly')
    setup(system) # make igblast executables and prepare dependency files 

def main():
    """
    finds files in the folder that has the file extension 
    """
    try:
        print('Looking for executables and required databases')
        for root,folder,files in os.walk(os.environ['HOME']):
            for file in files:
                if 'igblastn' in file and 'ncbi' in root and release in root:
                    igblastpath=root
                    igblastcmd=os.path.join(root,file)
                if cn + '_' + receptor + 'V.fasta' in file and set1 in root:
                    germpath=root
                if cn + '_V.phr' in file and 'internal_data' in root:
                    internal=root
                if file==cn + '_gl.aux' and 'optional_file' in root:
                    optional=root
        if igblastcmd and optional and internal and germpath: # previous download good to go
            print('Found\n%s\n%s\n%s\n%s' %(igblastcmd,optional,internal,germpath))
            df=create_fasta_file(path,outname)
            runigblast(igblastcmd,germpath,optional)
        else:# one or all required directories and executables is/are missing will perform flesh install optional or internal or germpath
            install()# run setup
    except Exception as e:
        print(e)

if __name__=='__main__':
    main()

#print('\nRunning igblastplus2.py %s/%s*dna*.fasta  %s/%s_igblast.log %s_%s' %(working_dir,outname,working_dir,outname,outname,project))  


#subprocess.call('igblastplus2.py %s/%s*dna*.fasta  %s/%s_igblast.log %s_%s' %(working_dir,outname,working_dir,outname,outname,project),shell=True)  

