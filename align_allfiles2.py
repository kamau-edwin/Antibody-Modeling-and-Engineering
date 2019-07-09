#! /usr/bin/env python
# original Written by Jules Jacobsen (jacobsen@ebi.ac.uk). Feel free to do whatever you like with this code.
# modified by Edwin kamau on 7.15.15
# Does alignment and grabs score terms from clustered pdbs after design or loop modeling
# outputs a csv file with score terms and rmsd that can then be used for score vs rmsd plots
# extensively modified by Robert L. Campbell (rlc1@queensu.ca)
#from __future__ import print_function
from pymol import cmd
import glob
import re
import operator
import sys
import matplotlib as mpl
mpl.use('AGG')
import matplotlib.pyplot as plt

  # Main pymol script outputs ensemble figures and csv file used to make plot 
def align_allfiles(native=None,files='*.pdb',name='n. ca',cutoff=2.0, cycles=5):
  """
  Aligns all models in a list of files to one target using align. 

  usage:
        where target specifies the model id you want to align all others against,
        and  cutoff and cycles are options passed to the align command.
        You can specify the files to load and align using a wildcard. You should specify same
	number of residues for the target and other models. if outfile is not given output is
        stddount

    	By default all residues in the structure are aligned,cutoff is 2,number of cycles is 5,        and align method is used.

    Example:
      pymol -cdq "run /$SCRIPTS/align_allfiles;align_allfiles <target name>,[none default 
      arguments eg res_sel= i. 100-120,...]"

  """
  help='pymol -cdq "run/$SCRIPTS/align_allfiles;align_allfiles native=path to native> name="n. ca",cutoff=2,cycles=5 method=align'
  print native
  if not native:
    print 'Native structure not given will use lowest energy model.'
  if not files:
    print 'Missing pdb files or zipped pdb file.See Usage'
    print ''
    print help
    cmd.quit()
  file_list=glob.glob('*.pdb')
  #print file_list
#Define rmsd and residue selection for mobile and target
  all='all' 
  ca='n. ca'
  bb='n. n+c+ca+o'
  id=file_list[0][:5]
  print id
#obtain files from folder and remove native
  main="";dic={};energy=[]
  extension = re.compile( '(^.*[\/]|\.(pdb|ent|brk))' )
  for file in file_list:
    with open(file,'r') as f:
     for line in f:
	if line.startswith('pose'):
	  score=float('%0.2f'%(float(line.split()[-1])))
          energy.append(int(score))
 	  dic[file]=score
	  main=sorted(dic.items(),key=operator.itemgetter(1))
  print main
  low_score=min(energy)-1
  high_score=max(energy)+1
  fig,ax=plt.subplots()
  for pdb,score in main:
   if native:
     target=native
     cmd.load(target)
   else:
    target=main[0][0]
    cmd.load(target)
   target=extension.sub('',target)
   obj_name=extension.sub('',pdb)
   cmd.load(pdb,obj_name)
   rms_bb = cmd.align('%s & %s'%(obj_name,bb),'%s & %s'%(target,bb),cutoff=cutoff,cycles=cycles)
   rms_all = cmd.align('%s & %s'%(obj_name,all),'%s & %s'%(target,all),cutoff=cutoff,cycles=cycles)
   if obj_name==extension.sub('',main[0][0]):
       ax.scatter('%0.2f'%(rms_bb[0]),score,marker="*",color='red',s=25)
   else:
       ax.scatter('%0.2f'%(rms_bb[0]),score,marker="o",color='black')
   ax.set_title(id,fontweight='bold',fontsize=14,loc='center')
   ax.tick_params(axis='both',top='off',right='off',labelsize=8)
   #ax.set_xlim([-0.2,5])
   ax.set_ylim([low_score,high_score])
   ax.text(0.98, 0.04,'N=%s'%len(main),
           verticalalignment='bottom', horizontalalignment='right',
          transform=ax.transAxes,
          color='grey', fontsize=8)
  fig.text(0.5,0.052,'Backbone RMSD',ha='center',va='center',fontweight='semibold')
  fig.text(0.06, 0.5, 'SCORE', ha='center', va='center', rotation='vertical',fontweight='semibold')
  fig.savefig(id+'_backbone',dpi=300,bbox_inches='tight')  
  
cmd.extend('align_allfiles',align_allfiles)
