#!/usr/bin/python
########################################################################
# Perform binding energy analysis interanalyzer.py
#
# Last modified: 06.03.2015 by Edwin Kamau
#
# Usage: python Interface Analyzer 
########################################################################
# import all the necessary functions required to make functional calls
import string, sys, re, os
from subprocess import call
from random import randint

# If no arguments were given, print a helpful message

if len(sys.argv)!=3:
    print 'Usage: python Interanalyzer.py <PDBID> <chain(s) to keep fixed>'
    sys.exit(0)


pwd = os.getcwd()

# name that will be passed to all files and folders obtained from 
pdir = sys.argv[1]  # name of the job obtained from item 2 of usage(see above)
fixed=sys.argv[2]
dir=pdir[:-5]
mod=pdir[-4:]
#create project directory std output/error and 
#if not os.path.exists('results'): # project directory
#       os.makedirs('results')

#if not os.path.exists('scored_files'): # project directory
 #      os.makedirs('scored_files')

# Create flag files
#print 'Creating flag files...'
#f = open('./' + pdir + '.flags', 'w') # write below to pdir.flags in jobname folder
#specific options for InterfaceAnalyzer
#f.write('-in:file:s ' + pdir + '.pdb\n')
#f.write('-packstat true\n') #this will actually run rosetta's packstat calculation (slow)
#f.write('-interface ' + fix_chain + '\n')
#f.write('-tracer_data_print false\n') #make a score file with all the important info instead of just printing to the terminal
#f.write('-out:file:score_only\n') #This will cause output of all of the info to a file called "pack_input_score.sc"
#f.write('-out:file:scorefile %s_analysis.sc\n'% pdir[:-5])# % pdir) # path to scorefile  
#f.write('-out:file:scorefile results/interface_%s.sc\n' % pdir) # path to scorefile  
#f.write('-pack_input true\n') #will relax the input interface residues (good if you are starting from raw pdb files as in this case)
#f.write('-pack_separated true\n') #will also pack the monomers to calculated dG bind.
#f.write('-add_regular_scores_to_scorefile true\n') #will run the rest of rosetta's score function on your complex using score12
###these are some tweeks that we have found helpful
#f.write('-atomic_burial_cutoff 0.01\n') #This is set to help rosetta identify buried polar atoms properly
#f.write('-sasa_calculator_probe_radius 1.4\n') #This is the default water probe radius for SASA calculations, sometimes lowering the radius helps rosetta more accurately find buried polar atoms 
#f.write('-pose_metrics::interface_cutoff 6.0\n') # this defines how far away a CBeta atom can be from the other chain to be considered an interface residue
#f.write('-pose_metrics::inter_neighbor_interface_cutoff 8.0\n') # if using multichain this defines how far apart nbr are detected btn fixed group and separate group  away a CBeta atom can be from the other chain to be considered an interface residue
#f.write('-inter_group_neighbors_cutoff 8.0\n')
###options to help rosetta pack the input interfaces
#f.write('-use_input_sc\n') # will include the input rotamer in packing
#f.write('-ex1\n') #expand rotamer library around chi1 rotamers
#f.write('-ex2\n') #expand rotamer library around chi2 rotamers
#f.write('-extrachi_cutoff 1\n') #this will build extra rotamers at all positions and not just core
#f.write('-overwrite\n')
#f.close()

#Create a shell script bash file. Bash file called by qsub 
f = open(dir + '/' + pdir + '.bash', 'w')
f.write('#!/bin/bash\n') 
f.write('#$ -S /bin/bash\n') #use this shell to interpret script
f.write('#$ -N ' + pdir + '\n') #name of the job
f.write('#$ -cwd\n') # run script in current working directoryi
f.write('#$ -pe openmpi 10-11\n')
f.write('#$ -M etk243@nyumc.org\n') #send status updates
f.write('#$ -j y \n') #output input and output into the same file
f.write('#$ -o %s/error.out\n'% dir) # output file
f.write('cd %s\n' % dir)
f.write('module load openmpi/gcc/64/1.4.5\n')
f.write('#mkdir -p relax\n')
f.write('#mpirun -np $NSLOTS relax.mpi.linuxgccrelease -database $ROSETTA_DATABASE -s %s.pdb  -in:file:fullatom -ignore_unrecognized_res -relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains -relax:ramp_constraints false -use_input_sc -ex1 -ex2 -nstruct 10 -out:path:pdb relax -out:path:score relax -no_optH false -flip_HNQ -overwrite>relax.log\n' % pdir)
f.write('#InterfaceAnalyzer.linuxgccrelease -database $ROSETTA_DATABASE -s relax/%s_*.pdb -packstat true -interface %s -tracer_data_print false -out:file:score_only -out:file:scorefile %s_analysis.sc -pack_input true -pack_separated true -add_regular_scores_to_scorefile true -atomic_burial_cutoff 0.01 -use_input_sc -ex1 -ex2 -inter_group_neighbors_cutoff 5.0 -out:nooutput>inter.log\n'% (pdir,fixed,mod))
f.write("awk '{print$38,$2,$6,$7,$9,$25,$33,$28,$13,$17}' *_analysis.sc|sed 's/_00..//g'|sed '/^ *$/d'|sed '2,${/description/d}'|sed 's/ /,/g'>%s.csv\n" % dir)
#f.write('InterfaceAnalyzer.linuxgccrelease @' + pdir + '.flags -database $ROSETTA_DATABASE &>interface.log\n')
f.close()

#Create the Launcher script to run sge
f = open( 'job.launch', 'a')
f.write('qsub %s/%s.bash\n' % (dir,pdir))
f.close()
call("chmod +x job.launch", shell=True)

print ''
print 'NEXT STEPS: Run the following command to launch the parallel jobs'
print '        ./job.launch'
print ''


