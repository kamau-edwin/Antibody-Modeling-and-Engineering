#! /usr/bin/env python
import operator
import re
import glob
from pymol import cmd
from subprocess import call
import time
import sys
import os

def loop_model_analysis(pdbs=None,align_file=None,name='n. ca',cycles=5,cutoff=2.0):
  """
  Calculates backbone and all atom rmsd for the modeled loop region,obtains model score and modeled loop structural characteristics
  SCRIPT methodology:
  inputs pdbs in the current directory,reads each pdb file filtered by the cluster it belongs to and extracts the score.
  Creates a sorted by score cluster dictionary.The lowest scoring decoy from each cluster is then used as the target for backbone and all atom rmsd calculations based on the modeled loop region only
  Generates an output file containing decoys with helical, strand or mixed conformation, with their score, rmsd, and modeled loop secondary structure 
  Calculates statistics for each cluster and overall loop modeling  
  Defaults:
	pymol alignment cutoff and cycles. 
  Needed
    At least one selection 
    Alignment file path
    Decoys
  Example
    pymol -cqd "run path to analysis file;analysis pdbs= *.pdb,cutoff,cycles,..."
    """    
  help='Usage:pymol -cqd "run analysis.py;analysis pdbs= path_to_pdbs or tar zipped pdb file,align_file=path_to_alignment_file, name=n. ca, cutoff=2.0, cycles=5"\nNote: Default alignment cutoff and cycles can be changed if necessary'
  print ''
  if not pdbs:
    print 'Missing pdb files or zipped pdb file.See Usage'
    print ''
    print help
    cmd.quit()
  if pdbs.split('.')[1]=='pdb' and os.path.exists('c.0.0.pdb'):
     files=glob.glob(pdbs)
  elif pdbs.split('.')[1]=='pdb' and not os.path.exists('c.0.0.pdb'):
     print 'tarred file exist.Using it instead of given pdbs argument'
     pdbs='*.tgz'
     print 'untarring file...'
     call('tar -zxf %s' %pdbs,shell=True)
     files=glob.glob('*.pdb')
  elif pdbs.split('.')[1]=='tgz' and not os.path.exists('c.0.0.pdb'):
     print 'untarring file ...'
     call('tar -zxf %s' %pdbs,shell=True)
     files=glob.glob('*.pdb')
  else:
     files=glob.glob('*.pdb')
  if not align_file:
    print 'Target-template alignment file missing.See Usage'
    print ''
    print help
    cmd.quit()
  ali_file=glob.glob(align_file) # alignment file with gaps
  extension = re.compile( '(^.*[\/]|\.(pdb|ent|brk))' )
    # dictionary to hold decoy and score
  #Create dictionary to sequester different clusters
  clust0_dic={};clust1_dic={};clust2_dic={};clust3_dic={};clust4_dic={}
  # sort by energy variables
  cluster0="";cluster1="";Cluster2="";cluster3="";cluster4="";energy=[]
  # create model sel1 and sel2 files from alignment file
  beg=[] # Modeled loop region beginning index
  end=[] # Modeled loop region end index
  # open the file and extract beginning and end of gaps and append
  with open(ali_file[0],'r') as f:
    for line in f:
      if '-' in line:
        for m in re.finditer('-+',str(line)):
          beg.append(m.start()+1)
          end.append(m.end())
  if len(beg)==2 and len(end)==2: # check if more than one loop was modeled
      sel1='i. %s-%s' %(beg[0],end[0])
      sel2='i. %s-%s' %(beg[1],end[1])
  elif len(beg)==1 and len(end)==1: # check if one loop was modeled and whether V1 or V2
    if beg[0]<34:
      sel1='i. %s-%s' %(beg[0],end[0])
      sel2=None
    elif beg[0]>34:
      sel2='i. %s-%s' %(beg[0],end[0])
      sel1=None
    # obtain modeled residues from modeling selection 
  print 'Selection 1: %s' %sel1
  print 'Selection 2: %s\n' %sel2
# extract score for each decoy create a dictionary and sort it
  for pdb in files:
    obj_name=extension.sub('',pdb)
    cmd.load(pdb,obj_name)    
    with open(pdb,'r') as f:
      for line in f:
        if line.startswith('pose'):
          score=float('%0.2f'%float(line.split()[-1]))
          energy.append(int(score))
          if pdb[2]=='0':
            clust0_dic[pdb]= score
            cluster0=sorted(clust0_dic.items(),key=operator.itemgetter(1))
          elif pdb[2]=='1':
            clust1_dic[pdb]=score
            cluster1=sorted(clust1_dic.items(),key=operator.itemgetter(1))
          elif pdb[2]=='2':
              clust2_dic[pdb]= score
              Cluster2=sorted(clust2_dic.items(),key=operator.itemgetter(1))

            cluster3=sorted(clust3_dic.items(),key=operator.itemgetter(1))
          elif pdb[2]=='4':
            clust4_dic[pdb]= score
            cluster4=sorted(clust4_dic.items(),key=operator.itemgetter(1))
          else:
            print 'No file with decoys found or clusters exceed set limit'
            print help
# RUN SECONDARY STRUCTURE ASSIGNMENT
  print 'Assigning Secondary structure...'
  if sel1 :
    print 'Processing based on selection 1'
    #log=open('cluster0.csv','a')
    v1_1H_c0=[];v1_2H_c0=[];v1_3H_c0=[];v1_4H_c0=[]
    v1_1S_c0=[];v1_2S_c0=[];v1_3S_c0=[];v1_4S_c0=[]
    v1_1HS_c0=[];v1_2HS_c0=[];v1_3HS_c0=[];v1_4HS_c0=[]
    for model,score in cluster0:
      obj_name=extension.sub('',model)
      cluster={'sec':[]}
      cmd.iterate('%s & %s & %s'%(obj_name,sel1,name),'sec.append(ss)',space=cluster)
      cluster=''.join(cluster['sec'])
      if 'H' in cluster and 'S' not in cluster:
        ss=[m.start() for m in re.finditer('H',str(cluster))]
        if len(ss)==1:
          v1_1H_c0.append(model)
        if len(ss)==2:
          v1_2H_c0.append(model)
        if len(ss)==3:
          v1_3H_c0.append(model)
        if len(ss)>=4:
          v1_4H_c0.append(model)
      if 'S' in cluster and 'H' not in cluster:
        ss=[m.start() for m in re.finditer('S',str(cluster))]
        if len(ss)==1:
            v1_1S_c0.append(model)
        if len(ss)==2:
            v1_2S_c0.append(model)
        if len(ss)==3:
            v1_3S_c0.append(model)
        if len(ss)>=4:
            v1_4S_c0.append(model)
      if 'H' in cluster and 'S' in cluster:
        ss=[m.start() for m in re.finditer('H|S',str(cluster))]
        if len(ss)==1:
            v1_1HS_c0.append(model)
        if len(ss)==2:
            v1_2HS_c0.append(model)
        if len(ss)==3:
            v1_3HS_c0.append(model)
        if len(ss)>=4:
            v1_4HS_c0.append(model)
# CLUSTER 1 Sel1 only
    v1_1H_c1=[];v1_2H_c1=[];v1_3H_c1=[];v1_4H_c1=[]
    v1_1S_c1=[];v1_2S_c1=[];v1_3S_c1=[];v1_4S_c1=[]
    v1_1HS_c1=[];v1_2HS_c1=[];v1_3HS_c1=[];v1_4HS_c1=[]
    for model,score in cluster1:
      obj_name=extension.sub('',model)
      cluster={'sec':[]}
      cmd.iterate('%s & %s & %s'%(obj_name,sel1,name),'sec.append(ss)',space=cluster)
      cluster=''.join(cluster['sec'])
      if 'H' in cluster and 'S' not in cluster:
        ss=[m.start() for m in re.finditer('H',str(cluster))]
        if len(ss)==1:
          v1_1H_c1.append(model)
        if len(ss)==2:
          v1_2H_c1.append(model)
        if len(ss)==3:
          v1_3H_c1.append(model)
        if len(ss)>=4:
          v1_4H_c1.append(model)
      if 'S' in cluster and 'H' not in cluster:
        ss=[m.start() for m in re.finditer('S',str(cluster))]
        if len(ss)==1:
          v1_1S_c1.append(model)
        if len(ss)==2:
          v1_2S_c1.append(model)
        if len(ss)==3:
          v1_3S_c1.append(model)
        if len(ss)>=4:
          v1_4S_c1.append(model)
      if 'H' in cluster and 'S' in cluster:
        ss=[m.start() for m in re.finditer('H|S',str(cluster))]
        if len(ss)==1:
          v1_1HS_c1.append(model)
        if len(ss)==2:
          v1_2HS_c1.append(model)
        if len(ss)==3:
          v1_3HS_c1.append(model)
        if len(ss)>=4:
          v1_4HS_c1.append(model)
# CLUSTER 2 Sel1 only
    v1_1H_c2=[];v1_2H_c2=[];v1_3H_c2=[];v1_4H_c2=[]
    v1_1S_c2=[];v1_2S_c2=[];v1_3S_c2=[];v1_4S_c2=[]
    v1_1HS_c2=[];v1_2HS_c2=[];v1_3HS_c2=[];v1_4HS_c2=[]
    for model,score in Cluster2:
      obj_name=extension.sub('',model)
      cluster={'sec':[]}
      cmd.iterate('%s & %s & %s'%(obj_name,sel1,name),'sec.append(ss)',space=cluster)
      cluster=''.join(cluster['sec'])
      if 'H' in cluster and 'S' not in cluster:
        ss=[m.start() for m in re.finditer('H',str(cluster))]
        if len(ss)==1:
          v1_1H_c2.append(model)
        if len(ss)==2:
          v1_2H_c2.append(model)
        if len(ss)==3:
          v1_3H_c2.append(model)
        if len(ss)>=4:
          v1_4H_c2.append(model)
      if 'S' in cluster and 'H' not in cluster:
        ss=[m.start() for m in re.finditer('S',str(cluster))]
        if len(ss)==1:
          v1_1S_c2.append(model)
        if len(ss)==2:
          v1_2S_c2.append(model)
        if len(ss)==3:
          v1_3S_c2.append(model)
        if len(ss)>=4:
          v1_4S_c2.append(model)
      if 'H' in cluster and 'S' in cluster:
        ss=[m.start() for m in re.finditer('H|S',str(cluster))]
        if len(ss)==1:
          v1_1HS_c2.append(model)
        if len(ss)==2:
          v1_2HS_c2.append(model)
        if len(ss)==3:
          v1_3HS_c2.append(model)
        if len(ss)>=4:
          v1_4HS_c2.append(model)
# CLUSTER 3 Sel1 only
    v1_1H_c3=[];v1_2H_c3=[];v1_3H_c3=[];v1_4H_c3=[]
    v1_1S_c3=[];v1_2S_c3=[];v1_3S_c3=[];v1_4S_c3=[]
    v1_1HS_c3=[];v1_2HS_c3=[];v1_3HS_c3=[];v1_4HS_c3=[]
    for model, score in cluster3:
      obj_name=extension.sub('',model)
      cluster={'sec':[]}
      cmd.iterate('%s & %s & %s'%(obj_name,sel1,name),'sec.append(ss)',space=cluster)
      cluster=''.join(cluster['sec'])
      if 'H' in cluster and 'S' not in cluster:
        ss=[m.start() for m in re.finditer('H',str(cluster))]
        if len(ss)==1:
          v1_1H_c3.append(model)
        if len(ss)==2:
          v1_2H_c3.append(model)
        if len(ss)==3:
          v1_3H_c3.append(model)
        if len(ss)>=4:
          v1_4H_c3.append(model)
      if 'S' in cluster and 'H' not in cluster:
        ss=[m.start() for m in re.finditer('S',str(cluster))]
        if len(ss)==1:
          v1_1S_c3.append(model)
        if len(ss)==2:
          v1_2S_c3.append(model)
        if len(ss)==3:
          v1_3S_c3.append(model)
        if len(ss)>=4:
          v1_4S_c3.append(model)
      if 'H' in cluster and 'S' in cluster:
        ss=[m.start() for m in re.finditer('H|S',str(cluster))]
        if len(ss)==1:
          v1_1HS_c3.append(model)
        if len(ss)==2:
          v1_2HS_c3.append(model)
        if len(ss)==3:
          v1_3HS_c3.append(model)
        if len(ss)>=4:
          v1_4HS_c3.append(model)
  # CLUSTER  4  Sel1 only
    v1_1H_c4=[];v1_2H_c4=[];v1_3H_c4=[];v1_4H_c4=[]
    v1_1S_c4=[];v1_2S_c4=[];v1_3S_c4=[];v1_4S_c4=[]
    v1_1HS_c4=[];v1_2HS_c4=[];v1_3HS_c4=[];v1_4HS_c4=[]
    for model,score in cluster4:
      obj_name=extension.sub('',model)
      cluster={'sec':[]}
      cmd.iterate('%s & %s & %s'%(obj_name,sel1,name),'sec.append(ss)',space=cluster)
      cluster=''.join(cluster['sec'])
      if 'H' in cluster and 'S' not in cluster:
        ss=[m.start() for m in re.finditer('H',str(cluster))]
        if len(ss)==1:
          v1_1H_c4.append(model)
        if len(ss)==2:
          v1_2H_c4.append(model)
        if len(ss)==3:
          v1_3H_c4.append(model)
        if len(ss)>=4:
          v1_4H_c4.append(model)
      if 'S' in cluster and 'H' not in cluster:
        ss=[m.start() for m in re.finditer('S',str(cluster))]
        if len(ss)==1:
          v1_1S_c4.append(model)
        if len(ss)==2:
          v1_2S_c4.append(model)
        if len(ss)==3:
          v1_3S_c4.append(model)
        if len(ss)>=4:
          v1_4S_c4.append(model)
      if 'H' in cluster and 'S' in cluster:
        ss=[m.start() for m in re.finditer('H|S',str(cluster))]
        if len(ss)==1:
          v1_1HS_c4.append(model)
        if len(ss)==2:
          v1_2HS_c4.append(model)
        if len(ss)==3:
          v1_3HS_c4.append(model)
        if len(ss)>=4:
          v1_4HS_c4.append(model)
    log=open('v1_bysecstruct.txt','w')
    input_list=[v1_1H_c0,v1_2H_c0,v1_3H_c0,v1_4H_c0,v1_1S_c0,v1_2S_c0,v1_3S_c0,v1_4S_c0,v1_1HS_c0,v1_2HS_c0,v1_3HS_c0,v1_4HS_c0,v1_1H_c1,v1_2H_c1,v1_3H_c1,v1_4H_c1,v1_1S_c1,v1_2S_c1,v1_3S_c1,v1_4S_c1,v1_1HS_c1,v1_2HS_c1,v1_3HS_c1,v1_4HS_c1,v1_1H_c2,v1_2H_c2,v1_3H_c2,v1_4H_c2,v1_1S_c2,v1_2S_c2,v1_3S_c2,v1_4S_c2,v1_1HS_c2,v1_2HS_c2,v1_3HS_c2,v1_4HS_c2,v1_1H_c3,v1_2H_c3,v1_3H_c3,v1_4H_c3,v1_1S_c3,v1_2S_c3,v1_3S_c3,v1_4S_c3,v1_1HS_c3,v1_2HS_c3,v1_3HS_c3,v1_4HS_c3,v1_1H_c4,v1_2H_c4,v1_3H_c4,v1_4H_c4,v1_1S_c4,v1_2S_c4,v1_3S_c4,v1_4S_c4,v1_1HS_c4,v1_2HS_c4,v1_3HS_c4,v1_4HS_c4]
    for l in input_list:
      if not l:continue
      print>>log,l[0]
    log.close()
 #SELECTION TWO ONLY 
  if sel2:
    print 'processing based on selection 2'
    v2_1H_c0=[];v2_2H_c0=[];v2_3H_c0=[];v2_4H_c0=[]
    v2_1S_c0=[];v2_2S_c0=[];v2_3S_c0=[];v2_4S_c0=[]
    v2_1HS_c0=[];v2_2HS_c0=[];v2_3HS_c0=[];v2_4HS_c0=[]
    for model,score in cluster0:
      obj_name=extension.sub('',model)
      cluster={'sec':[]}
      cmd.iterate('%s & %s & %s'%(obj_name,sel2,name),'sec.append(ss)',space=cluster)
      cluster=''.join(cluster['sec'])
      if 'H' in cluster and 'S' not in cluster:
        ss=[m.start() for m in re.finditer('H',str(cluster))]
        if len(ss)==1:
          v2_1H_c0.append(model)
        if len(ss)==2:
          v2_2H_c0.append(model)
        if len(ss)==3:
          v2_3H_c0.append(model)
        if len(ss)>=4:
          v2_4H_c0.append(model)
      if 'S' in cluster and 'H' not in cluster:
        ss=[m.start() for m in re.finditer('S',str(cluster))]
        if len(ss)==1:
            v2_1S_c0.append(model)
        if len(ss)==2:
            v2_2S_c0.append(model)
        if len(ss)==3:
            v2_3S_c0.append(model)
        if len(ss)>=4:
            v2_4S_c0.append(model)
      if 'H' in cluster and 'S' in cluster:
        ss=[m.start() for m in re.finditer('H|S',str(cluster))]
        if len(ss)==1:
            v2_1HS_c0.append(model)
        if len(ss)==2:
            v2_2HS_c0.append(model)
        if len(ss)==3:
            v2_3HS_c0.append(model)
        if len(ss)>=4:
            v2_4HS_c0.append(model)
# CLUSTER 1 Sel2 only
    v2_1H_c1=[];v2_2H_c1=[];v2_3H_c1=[];v2_4H_c1=[]
    v2_1S_c1=[];v2_2S_c1=[];v2_3S_c1=[];v2_4S_c1=[]
    v2_1HS_c1=[];v2_2HS_c1=[];v2_3HS_c1=[];v2_4HS_c1=[]
    for model,score in cluster1:
      obj_name=extension.sub('',model)
      cluster={'sec':[]}
      cmd.iterate('%s & %s & %s'%(obj_name,sel2,name),'sec.append(ss)',space=cluster)
      cluster=''.join(cluster['sec'])
      if 'H' in cluster and 'S' not in cluster:
        ss=[m.start() for m in re.finditer('H',str(cluster))]
        if len(ss)==1:
          v2_1H_c1.append(model)
        if len(ss)==2:
          v2_2H_c1.append(model)
        if len(ss)==3:
          v2_3H_c1.append(model)
        if len(ss)>=4:
          v2_4H_c1.append(model)
      if 'S' in cluster and 'H' not in cluster:
        ss=[m.start() for m in re.finditer('S',str(cluster))]
        if len(ss)==1:
          v2_1S_c1.append(model)
        if len(ss)==2:
          v2_2S_c1.append(model)
        if len(ss)==3:
          v2_3S_c1.append(model)
        if len(ss)>=4:
          v2_4S_c1.append(model)
      if 'H' in cluster and 'S' in cluster:
        ss=[m.start() for m in re.finditer('H|S',str(cluster))]
        if len(ss)==1:
          v2_1HS_c1.append(model)
        if len(ss)==2:
          v2_2HS_c1.append(model)
        if len(ss)==3:
          v2_3HS_c1.append(model)
        if len(ss)>=4:
          v2_4HS_c1.append(model)
# CLUSTER 2 Sel2 only
    v2_1H_c2=[];v2_2H_c2=[];v2_3H_c2=[];v2_4H_c2=[]
    v2_1S_c2=[];v2_2S_c2=[];v2_3S_c2=[];v2_4S_c2=[]
    v2_1HS_c2=[];v2_2HS_c2=[];v2_3HS_c2=[];v2_4HS_c2=[]
    for model,score in Cluster2:
      obj_name=extension.sub('',model)
      cluster={'sec':[]}
      cmd.iterate('%s & %s & %s'%(obj_name,sel2,name),'sec.append(ss)',space=cluster)
      cluster=''.join(cluster['sec'])
      if 'H' in cluster and 'S' not in cluster:
        ss=[m.start() for m in re.finditer('H',str(cluster))]
        if len(ss)==1:
          v2_1H_c2.append(model)
        if len(ss)==2:
          v2_2H_c2.append(model)
        if len(ss)==3:
          v2_3H_c2.append(model)
        if len(ss)>=4:
          v2_4H_c2.append(model)
      if 'S' in cluster and 'H' not in cluster:
        ss=[m.start() for m in re.finditer('S',str(cluster))]
        if len(ss)==1:
          v2_1S_c2.append(model)
        if len(ss)==2:
          v2_2S_c2.append(model)
        if len(ss)==3:
          v2_3S_c2.append(model)
        if len(ss)>=4:
          v2_4S_c2.append(model)
      if 'H' in cluster and 'S' in cluster:
        ss=[m.start() for m in re.finditer('H|S',str(cluster))]
        if len(ss)==1:
          v2_1HS_c2.append(model)
        if len(ss)==2:
          v2_2HS_c2.append(model)
        if len(ss)==3:
          v2_3HS_c2.append(model)
        if len(ss)>=4:
          v2_4HS_c2.append(model)
# CLUSTER 3 Sel2 only
    v2_1H_c3=[];v2_2H_c3=[];v2_3H_c3=[];v2_4H_c3=[]
    v2_1S_c3=[];v2_2S_c3=[];v2_3S_c3=[];v2_4S_c3=[]
    v2_1HS_c3=[];v2_2HS_c3=[];v2_3HS_c3=[];v2_4HS_c3=[]
    for model, score in cluster3:
      obj_name=extension.sub('',model)
      cluster={'sec':[]}
      cmd.iterate('%s & %s & %s'%(obj_name,sel2,name),'sec.append(ss)',space=cluster)
      cluster=''.join(cluster['sec'])
      if 'H' in cluster and 'S' not in cluster:
        ss=[m.start() for m in re.finditer('H',str(cluster))]
        if len(ss)==1:
          v2_1H_c3.append(model)
        if len(ss)==2:
          v2_2H_c3.append(model)
        if len(ss)==3:
          v2_3H_c3.append(model)
        if len(ss)>=4:
          v2_4H_c3.append(model)
      if 'S' in cluster and 'H' not in cluster:
        ss=[m.start() for m in re.finditer('S',str(cluster))]
        if len(ss)==1:
          v2_1S_c3.append(model)
        if len(ss)==2:
          v2_2S_c3.append(model)
        if len(ss)==3:
          v2_3S_c3.append(model)
        if len(ss)>=4:
          v2_4S_c3.append(model)
      if 'H' in cluster and 'S' in cluster:
        ss=[m.start() for m in re.finditer('H|S',str(cluster))]
        if len(ss)==1:
          v2_1HS_c3.append(model)
        if len(ss)==2:
          v2_2HS_c3.append(model)
        if len(ss)==3:
          v2_3HS_c3.append(model)
        if len(ss)>=4:
          v2_4HS_c3.append(model)
  # CLUSTER  4  Sel2 only
    v2_1H_c4=[];v2_2H_c4=[];v2_3H_c4=[];v2_4H_c4=[]
    v2_1S_c4=[];v2_2S_c4=[];v2_3S_c4=[];v2_4S_c4=[]
    v2_1HS_c4=[];v2_2HS_c4=[];v2_3HS_c4=[];v2_4HS_c4=[]
    for model,score in cluster4:
      obj_name=extension.sub('',model)
      cluster={'sec':[]}
      cmd.iterate('%s & %s & %s'%(obj_name,sel2,name),'sec.append(ss)',space=cluster)
      cluster=''.join(cluster['sec'])
      if 'H' in cluster and 'S' not in cluster:
        ss=[m.start() for m in re.finditer('H',str(cluster))]
        if len(ss)==1:
          v2_1H_c4.append(model)
        if len(ss)==2:
          v2_2H_c4.append(model)
        if len(ss)==3:
          v2_3H_c4.append(model)
        if len(ss)>=4:
          v2_4H_c4.append(model)
      if 'S' in cluster and 'H' not in cluster:
        ss=[m.start() for m in re.finditer('S',str(cluster))]
        if len(ss)==1:
          v2_1S_c4.append(model)
        if len(ss)==2:
          v2_2S_c4.append(model)
        if len(ss)==3:
          v2_3S_c4.append(model)
        if len(ss)>=4:
          v2_4S_c4.append(model)
      if 'H' in cluster and 'S' in cluster:
        ss=[m.start() for m in re.finditer('H|S',str(cluster))]
        if len(ss)==1:
          v2_1HS_c4.append(model)
        if len(ss)==2:
          v2_2HS_c4.append(model)
        if len(ss)==3:
          v2_3HS_c4.append(model)
        if len(ss)>=4:
          v2_4HS_c4.append(model)
    #Check secondary structure and output lowest structure
    log=open('v2_bysecstruct.txt','w')
    input_list2=[v2_1H_c0,v2_2H_c0,v2_3H_c0,v2_4H_c0,v2_1S_c0,v2_2S_c0,v2_3S_c0,v2_4S_c0,v2_1HS_c0,v2_2HS_c0,v2_3HS_c0,v2_4HS_c0,v2_1H_c1,v2_2H_c1,v2_3H_c1,v2_4H_c1,v2_1S_c1,v2_2S_c1,v2_3S_c1,v2_4S_c1,v2_1HS_c1,v2_2HS_c1,v2_3HS_c1,v2_4HS_c1,v2_1H_c2,v2_2H_c2,v2_3H_c2,v2_4H_c2,v2_1S_c2,v2_2S_c2,v2_3S_c2,v2_4S_c2,v2_1HS_c2,v2_2HS_c2,v2_3HS_c2,v2_4HS_c2,v2_1H_c3,v2_2H_c3,v2_3H_c3,v2_4H_c3,v2_1S_c3,v2_2S_c3,v2_3S_c3,v2_4S_c3,v2_1HS_c3,v2_2HS_c3,v2_3HS_c3,v2_4HS_c3,v2_1H_c4,v2_2H_c4,v2_3H_c4,v2_4H_c4,v2_1S_c4,v2_2S_c4,v2_3S_c4,v2_4S_c4,v2_1HS_c4,v2_2HS_c4,v2_3HS_c4,v2_4HS_c4]
    for l in input_list2:
      if not l:continue
      print>>log,l[0]
    log.close()
cmd.extend('analysis',loop_model_analysis)
