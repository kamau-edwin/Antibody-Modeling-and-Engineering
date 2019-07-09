#!/usr/bin/python
########################################################################
# File generator for comparative modeling protocol.py
#
# Written by: Edwin Kamau and Brett Spurrier 
# Last modified: 12.30.2014 by Edwin Kamau
#
# Usage: python Template.py project template  target_sequence template_sequence
########################################################################
# import all the necessary functions required to make functional calls
import string, sys, re, os
from subprocess import call
from random import randint

# def generate_flags_file(pdir, jobname, job):
    # '''
    # Given a project directory, job name, and job number, generates a flags file for Rosetta.
# 
    # '''
    
# If no arguments were given, print a helpful message
# Name assigned to fasta file, prefix for frag_files ensure it matches required flag calls  

if len(sys.argv)!=8:
    print 'Usage: python homology.py <project> <Template> <target_sequence> <template_sequence> <start> <end> <nstruct> <robetta job id>'
    sys.exit(0)
#global variables 
pwd = os.getcwd()
pdir = sys.argv[1]  # name of the job obtained from item 2 of usage(see above)
PDB = sys.argv[2]
start = int(sys.argv[5])
end = int(sys.argv[6])
nstruct=int(sys.argv[7])
#id = int(sys.argv[8])
#scripts = "/ifs/home/etk243/HIV/software"

for i in sys.argv[3:4]: # item 3 of usage (see heading above)
    try: 
        target_sequence=i
    except ValueError:
        print 'alignment <target_sequence>'
        sys.exit(0)
for i in sys.argv[4:5]: # item 4 of usage (see heading above)
    try: 
        template_sequence=i
    except ValueError:
        print 'alignment <template_sequence>'
        sys.exit(0)
    else:
        print 'Aligning %s with clustalw2' % (target_sequence)
        f = open('./temp.fasta', 'w') 
        f.write('>Target\n%s\n>Template1\n%s' % (target_sequence, template_sequence))
        f.close()

#Run the alignment
#return_code = call("clustalo -i temp.fasta -o temp.aln --outfmt=clu", shell=True)  # maybe True works better?
return_code = call("clustalw2 temp.fasta",shell=True)

#Parse the alignment file
file = open("./temp.aln", 'r')
full_target_sequence = ''
full_template_sequence = ''

while 1:
    line = file.readline() # while loop that reads temp.aln line by line 
    if not line:
        break # stop while loop if no line is found
    if line.find('Target')==0: # find object in position one
        reline = re.sub(' +',' ',line) # substitute any space with single space
        reline_split = reline.split( ) #split reline on the space
        full_target_sequence += reline_split[1] #append position 2 of the reline_split to file full target ... 
    if line.find('Template1')==0:
        reline = re.sub(' +',' ',line)
        reline_split = reline.split( )
        full_template_sequence += reline_split[1]        

#create project directory std output/error and 
if not os.path.exists(pdir): # project directory
       os.makedirs(pdir)
if not os.path.exists(pdir + '/' + pdir + '.out'): # Rosetta output file
       os.makedirs(pdir + '/'+ pdir + '.out')
if not os.path.exists(pdir + '/frags'): # Rosetta output file
       os.makedirs(pdir + '/frags')
if not os.path.exists(pdir + '/' + pdir + '.out/Combined'): # Rosetta output file
       os.makedirs(pdir + '/'+ pdir +'.out/Combined')

f = open(pdir + '/'+  pdir + '.ali', 'w') # will create an alignment file used in modelling 
f.write('CLUSTAL 2.1 Multiple Sequence Alignments\n%s  %s\n%s_relax.pdb  %s' % (pdir, full_target_sequence, PDB, full_template_sequence))
f.close()

# fasta file of the target sequence passed to the flag file 
f = open(pdir + '/' + pdir + '.fasta', 'w')
f.write('>'+ pdir + '\n%s\n\n' % (target_sequence))
f.close()


#print pdir + ' ' + full_target_sequence   
print PDB + '_relax.pdb\n' + full_template_sequence
print ''
print pdir + '\n' + full_target_sequence      
#Create SS file
#V1V2_start_C = full_target_sequence.find('C')+1   
#V1V2_end_C = full_target_sequence.rfind('C')+1
print ' '
#for m in re.finditer('C',str(full_target_sequence)):
#	print 'AtomPair CB %s CB   HARMONIC     0.01' % (m.start()+1)
#	print ' '
#	print 'AtomPair SG %s SG   HARMONIC 2.03  0.01' % (m.start()+1)
f = open(pdir + '/' + pdir + '.loops','w')
#for m in re.finditer('-+',str(full_target_sequence)):
	#f.write('LOOP %i %i 0 0 0\n' % (m.start()-1,m.end()))
for m in re.finditer('-+',str(full_template_sequence)):
	f.write('LOOP %i %i 0 0 0\n' % (m.start()-1,m.end()))
	print 'LOOP %i %i 0 0 0' % (m.start()-1,m.end())	
f.close()
#V1V2_end_C = full_template_sequence.rfind('C')+1
#V1V2_end_C = full_template_sequence.rfind('C')+1
#V1V2_end_C = full_template_sequence.rfind('C')+1
#V1_start_C = full_template_sequence.find('C', V1V2_start_C + 1)+1
#V2_end_C =  full_template_sequence.find('C', V1_start_C + 1)+1

#f = open(pdir + '/' + pdir + '.ss', 'w')
#f.write('%s %s\n' % (V1V2_start_C, V1V2_end_C))  # C1_start, C1_end
#f.write('%s %s\n' % (V1_start_C, V2_end_C))      # C2_start, C2_end ???
#f.close()

#create full atom constraint for full atom S atom
#f = open(pdir + '/' + pdir + '.fa', 'w')
#f.write('AtomPair SG %s SG %s HARMONIC 2.03 0.1\n' % (V1V2_start_C, V1V2_end_C))  # C1_start, C1_end
#f.write('%s %s\n' % (V1_start_C, V2_end_C))      # C2_start, C2_end ???
#f.close()

# create centroid constraint file for cysteine atom 
#f = open(pdir + '/' + pdir + '.cst', 'w')
#f.write('AtomPair CB %s CB %s HARMONIC 4.0 0.1\n' % (V1V2_start_C, V1V2_end_C))  # C1_start, C1_end
#f.close()
#Cleanup
return_code = call("rm -f temp.*", shell=True)

# Create flag files
print 'Creating flag files...'
f = open(pdir + '/' + pdir + '.flags', 'w') # write below to pdir.flags in jobname folder
f.write('-run:protocol threading\n')
f.write('-in:file:fasta ' + pdir + '.fasta\n') # target sequence fasta file
f.write('-in:file:alignment ' + pdir + '.ali\n')
f.write('-in:file:template_pdb ../starting_files/start_input/'+ PDB + '_relax.pdb\n')       
f.write('-in:file:fullatom\n')
f.write('-cm:aln_format general\n')
#f.write('-in:file:native ../'+ PDB + '.pdb\n')
f.write('-in:file:psipred_ss2 frags/' + pdir + '.psipred_ss2\n')
#f.write('-fix_disulf ' + pdir + '.ss\n')
#f.write('-build_disulf\n' )# flag for remodel protocol
#f.write('-detect_disulfide_before_relax true\n')#flag for abinitio protocol
#f.write('-constraints:cst_fa_file ' + pdir + '.fa\n')
#f.write('-constraints:cst_fa_weight 10\n')
#f.write('-constraints:cst_file ' + pdir + '.cst\n')
#f.write('-constraints:cst_weight 10\n')
#f.write('-frag3 ' + pdir + '.frag3\n')
#f.write('-frag9 ' + pdir + '.frag9\n')
f.write('-loops:loop_file ' + pdir + '.loops\n')
#f.write('-loops:frag_files ' + pdir + '.frag9 ' + pdir + '.frag3 none\n')
f.write('-loops:frag_files frags/' + pdir + '.200.9mers frags/' + pdir + '.200.3mers none\n')
f.write('-loops:frag_sizes 9 3 1\n')
f.write('-idealize_after_loop_close\n')
f.write('-loops:extended\n')
#f.write('-loops:build_initial\n')
f.write('-loops:remodel quick_ccd\n')
f.write('-loops:refine refine_ccd\n')
#f.write('#-loops:refine refine_kic\n')
f.write('-loops:relax fastrelax\n')
#f.write('#-random_grow_loops_by 4\n')
f.write('-loops:fast\n')
f.write('#-relax:thorough\n')
f.write('-out:nstruct %i\n' % nstruct)
f.write('-out:file:silent '+ pdir + '.out/' + pdir + '.silent.${SGE_TASK_ID} \n')
f.write('-out:file:scorefile '+ pdir + '.out/threading_score.fasc\n')
f.write('-out:file:fullatom\n')
f.write('-out:file:silent_struct_type binary\n')
f.write('-use_input_sc\n')
f.write('-ex1\n')
f.write('-ex2\n')
f.write('-overwrite\n')
f.write('-run:seed_offset ${SGE_TASK_ID}\n')
f.close()

#Create a shell script bash file. Bash file called by qsub 
f = open(pdir + '/' + pdir + '.bash', 'w')
f.write('#!/bin/bash\n') 
f.write('#$ -S /bin/bash\n') #use this shell to interpret script
f.write('#$ -t %i-%i\n'%(start,end))
f.write('#$ -N ' + pdir + '\n') #name of the job
f.write('#$ -cwd\n') # run script in current working directory
f.write('#$ -M etk243@nyumc.org\n') #send status updates
f.write('#$ -j y \n') #output input and output into the same file
#f.write('#$ -o ' + os.getcwd() + '/' + pdir +'/error.out\n') # output file
f.write('#$ -o error.out\n') # output file
f.write('module load openmpi/gcc/64/1.4.5\n')
#f.write('cd ' + os.getcwd() + '/' + pdir + '\n') # output file
f.write('minirosetta.linuxgccrelease @' + pdir + '.flags -database $ROSETTA_DATABASE &>log.out\n')
f.write('exit 0;\n')
f.close()

#Create a clustering and score file 
f=open(pdir + '/' + pdir + '.out/Combined/post_processing.bash', 'w')
f.write('#!/bin/bash\n')
f.write('#$ -S /bin/bash\n')
f.write('#$ -N '+ pdir +'_anal\n')
f.write('#$ -cwd\n')
f.write('#$ -j y\n')
f.write('#$ -o '+ pdir + '/'+ pdir +'.out/Combined/processing.out\n')
f.write('module load openmpi/gcc/64/1.4.5\n')
f.write('cd '+ pdir + '/' + pdir +'.out/Combined\n')
f.write('combine_silent.linuxgccrelease -in:file:silent ../*.silent.* -out:file:silent Combined.sile -out:file:silent_struct_type binary -silent_read_through_errors -database $ROSETTA_DATABASE >&Combine.log\n')
f.write('Clusteranscore.py Combined.sile per.silent -1 0.1\n')
f.write('rm -f ../*.silent.*\n')
f.write('mkdir -p clust\n')
f.write('cd clust\n')
f.write('cluster.linuxgccrelease -in:file:silent ../*.silent -in:file:fullatom -in:file:native $INPUT/'+ PDB + '_relax.pdb -cluster:gdtmm -cluster:radius -1 -cluster:population_weight 0.0 -cluster:sort_groups_by_energy -no_optH -out:pdb -database $ROSETTA_DATABASE >& cluster.log\n')
f.write('mpirun -np 2 score.linuxgccrelease -in:file:s *.pdb -in:file:native $INPUT/'+ PDB + '_relax.pdb -out:file:scorefile score.sc -database $ROSETTA_DATABASE >score.log\n')

# move top10 structures based on score
f.write("#awk '{print$2,$30}' score.sc |sort -nk1 | awk '{print$2}' >byscore\n")
f.write("awk '{print$2,$21}' score.sc |sort -nk1 | awk '{print$2}' >byscore\n")
f.write("sed -i 's/_0001/.pdb/' byscore\n")
f.write('i=0\n')
f.write('mkdir -p best20\n')
f.write('while read line\n')
f.write('do\n')
f.write('    if [ "$i" -lt 20 ]; then\n')
f.write('        mv $line best20\n')
f.write('	((i++))\n')
f.write('    fi\n')
f.write('done < byscore\n')
f.write('exit 0;\n')
f.close()
call('chmod +x '+ pdir + '/' + pdir + '.out'+'/Combined/post_processing.bash',shell=True)

#Create the Launcher script to run sge
f = open(pdir + '/job.launch', 'w')
#f.write('qsub ' + pdir + '/' + pdir + '.bash\n')
f.write('qsub ' + pdir + '.bash\n')
f.close()
call('chmod +x '+ pdir + '/job.launch', shell=True)

# Create analysis launch file:
f = open( 'analysis.launch', 'a')
f.write('qsub '+ pdir + '/' + pdir + '.out/Combined/post_processing.bash\n')
f.close()
call("chmod +x analysis.launch", shell=True)

# Download fragment files from robetta
f =open('frag_picker.launch','a')
f.write('qsub '+ pdir + '/' + 'frags/frags.bash\n')
f.close()
call('chmod +x frag_picker.launch',shell=True)

f = open(pdir + '/frags/frags.bash','w')
f.write('#!/bin/bash\n')
f.write('#$ -S /bin/bash\n')
f.write('#$ -N '+ pdir +'_frags\n')
f.write('#$ -cwd\n')
f.write('#$ -j y\n')
f.write('#$ -o log.out\n')
f.write('cd ' + pdir + '/frags\n')
f.write('module load openmpi/gcc/64/1.4.5\n')
f.write('make_fragments.pl -id %s -frag_sizes 9,3 ../%s.fasta\n' %(pdir,pdir))
f.write('ls *.* |egrep -v  "*mers|*ss2|*.bash" |xargs rm\n') 
f.write('cd ..\n')
f.write('./job.launch\n')
f.write('rm -f job.launch\n')
f.write('exit;\n')
#f.write('curl http://www.robetta.org/downloads/fragments/%i/aat000_09_05.200_v1_3 >%s/%s.frag9\n' % (id,pdir,pdir))
#f.write('curl http://www.robetta.org/downloads/fragments/%i/aat000_03_05.200_v1_3 >%s/%s.frag3\n' % (id,pdir,pdir))
#f.write('curl http://www.robetta.org/downloads/fragments/%i/t000_.psipred_ss2 >%s/%s.psipred_ss2\n' % (id,pdir,pdir))
#f.close()
#call("chmod +x downloadfromrobetta.launch", shell=True)
print ''
print 'NEXT STEPS: Run the following command to launch the parallel jobs'
print '        ./frag_picker.launch'
print ''
