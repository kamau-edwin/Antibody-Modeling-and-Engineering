#! /local/apps/python/3.4.3/bin/python
#############################################################
# File to analyze cartesian ddg and residue_energy_breakdown#
# Written by: Edwin Kamau                                   # 
# created on: 12.07.2017                                    #
# updated on: 01.13.18                                      # 
# Usage: ddg_analysis.py  FUN project name                  #
#############################################################

# import all necessary functions

#import glob
import operator
import matplotlib as mpl
mpl.use('AGG')
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import sys
import os 
import csv
import re
import random 
#import optparse import OptionParser


longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}


# take in input 
if len(sys.argv) < 5:
    print ('Usage: ddg_analysis.py <function_to_run> <project_name> <folder> <file_extension>\n folder:location where the files with the given extension are located\nfunction to run:\n dg_score:computes effect of mutation on stability,plots ddg score for each mutant\n\n pairwise_energies:computes residue-residue pairwise energies for a given mutant\n\n ddg_bind:computes binding energy metrics requires interface analyzer protocol to have been run before and data stored in a csv file\nNOTE:\nrequires chothia numbered pdb that is the same length as the renumbered pdb\nused for running the protocol.\nThe order of the chains should be HL:A where A is antigen')
    sys.exit(0)
#for i in range(len(sys.argv)):
 #  print(sys.argv[i])
if len(sys.argv)>5:
  project = sys.argv[2]
  folder=sys.argv[3]
  file_extension=sys.argv[4]
  file_extension2=sys.argv[5]
else:
  project = sys.argv[2]
  folder=sys.argv[3]
  file_extension=sys.argv[4]

def snugdock_stats(project,analysis,file_extension,file_extension2):
	'''
	searches for files in the analysis folder and groups them by snugfiles and abangle files and creates names that will be used to load the files as csv and group them by strain
	# creates the path which to store output analysis folder
	#files: files to be imported, 
	#strain: grouping list based on peptide docked 
	'''

	from scipy import stats 
	import matplotlib as mpl
	mpl.use('agg')
	import matplotlib.pyplot as plt 
	import seaborn as sns
	import pandas as pd
	import numpy as np
	import os
	import sys 
	plt.style.use('ggplot')

	# 
	snug_files=[];ab_files=[];path='';strain=[]
	for root,dirs,files in os.walk(os.getcwd()):
		for file in files:
			if folder in root and file.startswith(file_extension):
			   peptide=file.partition('_vh')[0].split('-')[-1]
			   snug_files.append(os.path.join(root,file))
			   path=root
			   strain.append(peptide)
			if 'analysis' in root and file.startswith(file_extension2):
				ab_files.append(os.path.join(root,file))

	# Take both snug and abangle files and merged them on description and create global variable
	cols=['description','total_score', 'rms', 'Fnat', 'I_sc', 'Irms','conf_num1', 'conf_num2'] # columns of interest
	header=['description', 'HL', 'HC1', 'HC2', 'LC1', 'LC2', 'dc']
	gnames=[]
	fig,ax=plt.subplots(figsize=(10,8))
	if snug_files and ab_files:
		for rec in snug_files:
			for rec2 in ab_files:
				for ind in range(len(strain)):
					if strain[ind] in rec: 
						if strain[ind] in rec2:
							name=''.join('m_'+ os.path.split(rec2)[1].partition('.sc')[0].partition('_')[2].replace('-','_'))
							gnames.append(name)
							m=pd.read_csv(rec,skiprows=1,delim_whitespace=True,usecols=cols)
							if ind==0:
								high=m[(m.Fnat>=0.5) & (m.rms<=1.0) | (m.Irms<=1.0)]
								medium=m[(m.Fnat>=0.3) & (m.rms>1.0) & (m.rms<=5.0) & (m.Irms>1.0) & (m.Irms<5)]
								acceptable=m[(m.Fnat>=0.1) & (m.Fnat<=0.3) & (m.rms>5.0) & (m.rms<=10.0) & (m.Irms>2) & (m.Irms<5)]
								incorrect=m[(m.Fnat<0.1) | (m.rms>10.0) & (m.Irms>4.0)]
								ax.scatter(high.rms,high.I_sc,marker="o",color='red',label='high')
								ax.scatter(medium.rms,medium.I_sc,marker="h",color='green',label='medium')
								ax.scatter(acceptable.rms,acceptable.I_sc,marker="d",color='blue',label='acceptable')
								ax.scatter(incorrect.rms,incorrect.I_sc,marker="p",color='black',label='incorrect')
								handle,label=ax.get_legend_handles_labels()
								plt.xlabel('Antigen rmsd',fontsize=14,fontweight='semibold')
								plt.ylabel('Interface score',fontsize=14,fontweight='semibold')
								lg=plt.legend(handle,label,loc=2,scatterpoints=1,fontsize=14,title='model quality',bbox_to_anchor=(0.99,1.01))
								lg_ttle=lg.get_title()
								lg_ttle.set_fontsize(14)
								lg_ttle.set_fontweight('bold')
								fig.savefig(path+'/' + name + '_model_quality',bbox_inches='tight',dpi=300)
							m=m[(m.rms<2) & (m.Fnat>0.6)].sort_values('I_sc')
							m['mutant']=name
							m2=pd.read_csv(rec2,names=header,delim_whitespace=True)
							globals()[name]=pd.merge(m2,m,on='description')

		# Perform tstats with the best models based on rms and 
		f=open(path + '/'+ project + '_ttest_stats.csv','w')
		f.write('mut1,mut2,angle,p_value\n')
		abangles=['HL','HC1','HC2','LC1','LC2','dc']
		for i in range(0,3):
			for j in range(1,4):
				for key in abangles:
					if i==j==1 or i==j==2 or i==2 and j==1:continue
					t,p=stats.ttest_ind(globals()[gnames[i]][key],globals()[gnames[j]][key],equal_var=False)
					f.write('%s,%s,%s,%s\n'%(gnames[i],gnames[j],key,p))
					if p<0.05:
						f.write('\n')
						f.write('%s,%s,%s,%s\n'%(gnames[i],gnames[j],key,p))
		f.close()
		for dframe in gnames:
			corr=globals()[dframe].iloc[:,1:7].corr()
			f,ax=plt.subplots(figsize=(8,6))
			sns.heatmap(corr, mask=np.zeros_like(corr, dtype=np.bool), cmap=sns.diverging_palette(220, 10, as_cmap=True),
	            square=True, ax=ax)
			f.savefig(''.join(path + '/' + dframe + '_heatmap'),bbox_inches='tight',dpi=300)
	else:
		print("looks like there are no files to process\nExiting ...")
		sys.exit(0)

def score_plot(scores,mutations,ylabel,xlabel,output_name,title=None,labels=None):
    """
    plots the results of the scoring terms dg sorted low to high.
    as computed by dg_score function.Assumes ddg has been calculated 
    and the individual terms have been decomposed. 
    outputs plots for each mutation
    ddg : calculated delta delta G based on cartesian ddg output total score avg(total_mutant)-avg(total_wt) 
    dg_scores: Decomposed, averaged and sorted difference in scoring terms with contribution greater/less than 0.1 kcal/mol. 
    """
    
    #ddg=f'{ddg:.2f}' # comment if using python 2.7
    #ddg='%.2f' % ddg # uncomment if using python 2.7 or below
    if labels:
       colors=[]
       for label in range(len(labels)):
         r = lambda: random.randint(0,255)
         colors.append('#%02X%02X%02X' % (r(),r(),r()))
    else:
       colors='#999999' 
    fig=plt.figure(figsize=(8,6))
    p=plt.bar(range(len(mutations)),scores,color=colors,align='center')
    if title:
       plt.title(title,fontsize=18,fontweight='bold')
       plt.xticks(range(len(mutations)),mutations,fontweight='semibold',fontsize=14)
    else:
       plt.xticks(range(len(mutations)),mutations,fontweight='semibold',fontsize=10)
    plt.ylabel(ylabel,fontsize=14,fontweight='bold')
    plt.xlabel(xlabel,fontsize=12,fontweight='semibold')
    if labels:
       lg=plt.legend(p,labels,bbox_to_anchor=(1.01, 1.0),loc=2, title='Mutation',borderaxespad=0.)
       lg_title=lg.get_title()
       lg_title.set_fontsize(14)
       lg_title.set_fontweight('bold')
    #fig.text(0.5,0.82,mutation,fontsize=14,fontweight='bold')
    #fig.text(0.76,0.14,'ddg=%s'%ddg,fontsize=10,fontweight='bold')
    #plt.show()
    fig.savefig(output_name,dpi=300,bbox_inches='tight')
    print ('Done plotting for %s for files in %s folder'%(sys.argv[1], folder))

def ddg_bind(path,project,file_list,resid,chain_h,chain_l):
    '''
    obtain the effect to binding energy for mutation not directly involved in binding an epitope or ligan 
    computes the ddg of binding between the wild type and the mutant
    requires project name and the file name to end with ana.sc
    returns ddg plot and summary of ddg binding
    '''
    print('Running %s for files in the %s folder ...\n' %(sys.argv[1],folder)) 
    mut_scan={}
    for file in file_list:
      mutation=os.path.split(file)[1].split('-')[1]
      wt_res,pose,mut_res=re.split('(\d+)',mutation)
      pose=int(pose)
      pdbid=resid[pose]
      dat=open(file,'r').readlines()[1:]
      wt={};mut={}
      for line in dat:
          if 'WT' in line:
              wt[line.split(',')[0]]=[float(line.split(',')[1]) for line in dat[1:] if 'WT' in line]
          if mutation in line:
              mut[mutation]=[float(line.split(',')[1]) for line in dat[1:] if mutation in line]
      ddg=sum(mut[mutation])/len(mut[mutation])-sum(wt['WT'])/len(wt['WT'])
      ddg=float('%.3f'%(ddg))
      mut_scan[mutation]=(pose,pdbid,ddg)
    mut_scan=sorted(mut_scan.items(),key=operator.itemgetter(1))
    f=open(path + '/' + project + '_interface_ddg_score_summary.csv','w')
    print('Writing binding energy results')
    print('poseid\tpdbid\tmutation  ddg')
    f.write('poseid,pdbid,mutation,ddg\n')
    for score in mut_scan:
          mutate=score[0]
          pose_id=score[1][0]
          pdb_id=score[1][1]
          score=score[1][2]
          f.write('%s,%s,%s,%s\n' %(pose_id,pdb_id,mutate,score))
          print('%s \t%s \t%s\t %s' %(pose_id,pdb_id,mutate,score))
    f.close()
    print('Generating plots\n')
    mutants=[mut[1][1] for mut in mut_scan]
    scores=[score[1][2] for score in mut_scan]
    labels=[''.join(re.split('(\d+)',mut[0])[0] +re.split('(\d+)',mut[1][1])[1]+ re.split('(\d+)',mut[0])[2]) for mut in mut_scan]
    if folder.startswith('interface'):
         output_name=path +'/'+ project + '_revert_to_germ_interface_analysis' 
   #      if not os.path.isfile(output_name +'.png'):
         title='Effect of mutation to %s binding affinity' %s
         xlabel='Residue ID'
         ylabel='r$\Delta \Delta G_{bind}(Kcal/mole)$'
         output_name=path +'/'+ project + '_revert_to_germ_interface_analysis' 
         score_plot(scores,mutants,ylabel,xlabel,output_name,title,labels)
         #elif os.path.isfile(output_name + '.png'):
         # print('Looks like you have previously plotted your data.\nIf you have new or updated data you will have to either delete the file\n\
   #rm -f %s.png file\n      OR\n Resubmit the script with a different project name\n' %(output_name))
    elif folder.startswith('alascan_int'): 
       output_name=path +'/'+ project + '_alascan_interface_analysis' 
       title='Effect of Alanine Scanning Mutagenesis %s binding' %project
       xlabel='Residue ID'
       ylabel='r$\Delta \Delta G_{bind}(Kcal/mole)$'
       score_plot(scores,mutants,ylabel,xlabel,output_name,title,labels)
    
def pairwise_energies(path,project,file_list,resid,chain_h,chain_l):
    """
    Obtains residue_residue pairwise interaction for designated twobody energies. Requires a non rosetta renumbered pdb in the path i.e project.pdb. 
    Output pairwise interaction between the wt and mutated residue with neighboring residues that make +/- 0.5 kcal interaction contribution.
    score
    """
    print('Running %s for files in the %s folder ...\n' %(sys.argv[1],folder)) 
    # get header description from the first file
    header=open(file_list[0],'r').readline().split()
    #define score terms and their indexes 
    score_terms=['fa_atr','fa_elec','fa_sol','fa_rep','hbond_sc','hbond_sr_bb','hbond_lr_bb','hbond_bb_sc','lk_ball','lk_ball_iso']
      
    score_terms_index=[index for index,term in enumerate(header) if term in score_terms]
    #process file list 
    for file in file_list:
    #    print ('\nprocessing file=%s\n' %file)
        mutation=os.path.split(file)[1].split('-')[0]
        pose=int(re.split('(\d+)',mutation)[1]) # pose id
        pdbid=resid[pose]
        #print(mutation,pose,pdbid) 
        record=open(file,'r').readlines() # breakdown .out file reading 
        pose_ids=[] # store buffer for pairwise energies based on pose id been analysed
        res_pair_term=[] # pairwise terms that are significant 
        for index in score_terms_index:
            stabilizing=sorted(list(set([int(line.split()[5]) for line in record if 'onebody' not in line and line.split()[2]==str(pose) and float(line.split()[index])<=float(-0.5)])))
            stab2=[(int(line.split()[5]),longer_names[line.split()[7]]) for line in record if 'onebody' not in line and line.split()[2]==str(pose) and float(line.split()[index])<=float(-0.5)]
            destabilizing=sorted(list(set([int(line.split()[5]) for line in record if 'onebody' not in line and line.split()[2]==str(pose) and float(line.split()[index])>=float(0.5)])))
            if stabilizing:
                pose_ids.append(stabilizing)
                print(stabilizing)
                res_pair_term.append(header[index])
            if destabilizing:
                pose_ids.append(destabilizing)
                res_pair_term.append(header[index])
        #print(pose_ids)
        pose_ids=sorted(list(set([pos for poze in pose_ids for pos in poze])))
        res_pair_term=list(set(res_pair_term))
        # open output file and write
        f=open(path + '/'+ mutation + '_energy_breakdown_by_score_terms.csv','a')
        f.write('\n')
        f.write('poseid,pdbid,wt_score,mut_score\n')
    #    print('poseid,pdbid,wt_score,mut_score\n')
        for term in res_pair_term:
            f.write('score term: %s\n' % term)
         #   print('term: %s \n' %term)
            for index in score_terms_index:
                for pos in range(len(pose_ids)):
                    pdb_id=resid[pose_ids[pos]]
                    if header[index]==term:
                        wt=[float(line.split()[index]) for line in record[1:] if 'onebody' not in line and line.split()[5]==str(pose_ids[pos]) and line.split()[2]==str(pose) and 'WT' in line]
                        mt=[float(line.split()[index]) for line in record[1:] if 'onebody' not in line and line.split()[5]==str(pose_ids[pos]) and line.split()[2]==str(pose) and 'MUT' in line]
                        if sum(wt)==0.0 or sum(mt)==0.0:continue
                        wt='{0:.3f}'.format(sum(wt)/len(wt))
                        mt='{0:.3f}'.format(sum(mt)/len(mt))

                        f.write('%s,%s,%s,%s\n' %(pose_ids[pos],pdb_id,wt,mt))
#                        print('\n%s,%s,%s,%s\n' %(pose_ids[pos],pdb_id,wt,mt))
        f.close()
def dg_score(path,project,file_list,resid,chain_h,chain_l):
    """
    obtain ddg and decompose all energy terms for wild type and mutant and compute the average and ddg between them requires ddg file and project name 
    returns plot and csv file of ddg and dg of the score terms
    """
    print('Running %s for files in the %s folder ...\n' %(sys.argv[1],folder)) 
    # define scoring terms to restrict to 
    twobody_terms=['fa_atr','fa_rep','fa_elec','fa_sol','lk_ball','lk_ball_iso','hbond_sr_bb','hbond_lr_bb','hbond_bb_sc','hbond_sc','rama_prepro','cart_bonded','ref']
    stabilizing = {}
    destabilizing= {}
    all_mutants= {}
    for file in file_list:
        mutation=os.path.split(file)[1].split('.')[0]
        pose=int(re.split('(\d+)',mutation)[1])
        pdbid=resid[pose]
        record=open(file,'r').readlines()
        wt_scores={};mut_scores={};mut=''; dg_scores={};ddg=''
        for rec in record:
            for i in range(3,len(rec.split()),2):
                if 'MUT' in rec:
                  mut_name =rec.split()[2]
                  mut_scores[rec.split()[i-1]]=[float(rec.split()[i]) for rec in record if 'MUT_' in rec]
                if 'WT:' in rec:
                  wt_scores[rec.split()[i-1]]=[float(rec.split()[i]) for rec in record if 'WT:' in rec] 
        wt_scores['ddg']=wt_scores.pop('WT:')
        #wt_scores.pop('WT:',{'ddg':wt_scores['WT:']})
        mut_scores['ddg']=mut_scores.pop(mut_name)
        # compute dg and ddg 
        for key in wt_scores.keys():
          mut=sum(mut_scores[key])/len(mut_scores[key])
          wt=sum(wt_scores[key])/len(wt_scores[key])
          dg=mut-wt
          if key == 'ddg':
            ddg=float('{0:.3f}'.format(dg))
            #ddg=float(f'{dg:.3f}')
            #ddg=float('%0.3f' % dg)
          elif key.split(':')[0] in twobody_terms:# filter dg score terms based on twobody contributing terms and scores that are greater or less than 0.1 kcal/mol
            if float(dg) >=0.5 or float(dg) < -0.5:
                dg_scores[key.split(':')[0]]=float('{0:.3f}'.format(dg))
                #dg_scores[key.split(':')[0]]=float('%0.3f'% dg) python 2.7 and below
        print('%s : ddg for %s mutation is %s' %(project,mutation,ddg)) 
        dg_scores=sorted(dg_scores.items(),key=operator.itemgetter(1))
        scores=[score[1] for score in dg_scores]
        terms=[term[0] for term in dg_scores]
        output_name=path + '/'+ project+'_'+ mutation+'_scoreterms'
        ylabel='$\Delta Energy (Kcal/mole)$'
        xlabel='score term'
        text='%s:%s' % (project,pdbid)
        score_plot(scores,terms,ylabel,xlabel,output_name)
        # write dg scores into a file as csv format
        #print (dg_scores)
        f=open(path + '/'+ mutation + '_energy_breakdown_by_score_terms.csv','w')
        outfile = csv.writer(f)
        for i in range(len(dg_scores[0])):
            outfile.writerow([x[i] for x in dg_scores])
        f.close()
        # get the two highest contributing terms to ddg
        if ddg > 0:
           destabilizing[mutation]=(int(pose),pdbid,ddg)
           all_mutants[mutation]=(int(pose),pdbid,ddg)
        elif ddg < 0:
           stabilizing[mutation]=(int(pose),pdbid,ddg)
           all_mutants[mutation]=(int(pose),pdbid,ddg)
    destabilizing=sorted(destabilizing.items(),key=operator.itemgetter(1))
    stabilizing=sorted(stabilizing.items(),key=operator.itemgetter(1))
    all_mutants=sorted(all_mutants.items(),key=operator.itemgetter(1))
    # 
    ddgs=[score[1][2] for score in all_mutants]
    spm=[mutant[1][1] for mutant in all_mutants]
    labels=[''.join(re.split('(\d+)',mut[0])[0] +re.split('(\d+)',mut[1][1])[1]+ re.split('(\d+)',mut[0])[2]) for mut in all_mutants]
    #LABELS
    output_name=path + '/'+ project+'_'+folder+'_dg_summary'
    #if not os.path.isfile(output_name+'.png'):
    ylabel='$\Delta \Delta G(Kcal/mole)$'
    title='Effect of mutation on antibody %s stability' % project
    xlabel='Residue ID' #(chothia numbering scheme)'
    score_plot(ddgs,spm,ylabel,xlabel,output_name,title,labels)
       
    # MAKE SUMMARY FILE FOR ALL MODELS
    f_all=open(path +'/'+ project + '_dg_summary.csv','w')
    f_all.write('Stabilizing mutants\n')
    f_all.write('Mutant,pose_id,pdb_id,ddg\n')
    for mut in stabilizing:
       f_all.write('%s,%s,%s,%s\n'%(mut[0],mut[1][0],mut[1][1],mut[1][2]))
    f_all.write('\n')
    f_all.write('Destabilizing mutants\n')
    f_all.write('Mutant,pose_id,pdb_id,ddg\n')
    for mt in destabilizing:
       f_all.write('%s,%s,%s,%s\n'%(mt[0],mt[1][0],mt[1][1],mt[1][2]))
    f_all.close()

def file_util(project,folder,file_extension):
    """
    finds files in the folder that has the file extension 
    returns, the path,project,file_list,len of chain 1 2 and 3
    (if available) 
    """
    print("entered file_util with %s %s %s" %(project,folder,file_extension))
    file_list=[]
    path=''
    pdb=''
    for root,dirs,files in os.walk(os.getcwd()):
        for file in files:
            if folder in root and file.endswith(file_extension):
                file_list.append(os.path.join(root,file))
                path=root
            if file ==''.join(project+'_chothia.pdb'):
                pdb=os.path.join(root,file)
    # get the first chain from pdb
    try:
        resid={}
        index=1
        pdb_ = open(pdb,'r').readlines()
        for line in pdb_:
            if 'CA' in line:
              chain=line.split()[4]
              res=line.split()[5]
              resid[index]=''.join(chain+res)
              index +=1
        chain_h=[resid[ind] for ind in range(1,len(resid)) if resid[ind].startswith('H')]
        chain_l=[resid[ind] for ind in range(1,len(resid)) if resid[ind].startswith('L')]
   #     for tup in range(1,len(resid)):
    #        print(tup,resid[tup])
        #print(chain_h,chain_l)
    except FileNotFoundError:
        print('pdb file not found in the path or its chain less \n')
        print('Exiting... please provide a file to assign pdbid')
        sys.exit(0)
    globals()[sys.argv[1]](path,project,file_list,resid,chain_h,chain_l)
if sys.argv[1] !='snugdock_stats':
   file_util(project,folder,file_extension)
if sys.argv[1] =='snugdock_stats':
   snugdock_stats(project,folder,file_extension,file_extension2)

#if __name__ == '--main__':
#Z	globals()[sys.argv[1]](project)

#Incase you run into unequel length of H chain and light chain has missing residues
#if pose > len(chain_h) and int(chain_h[-1].partition('H')[2]) > 112 and chain_l[0] !='L1':
 #          reduceby=int(chain_h[-1].partition('H')[2])-112# shift L residue by
  #        pdbid=resid[pose-reduceby] #pdb id
   #        print('chain H end residue is %s and therefore longer than usual fv length of 112\nAll chain L residues pose number will be reduced by %s to rectify the additional residue discrepancy with chothia numbering of chain L' % (chain_h[-1],reduceby))
