#!/usr/bin/python
########################################################################
# Perform binding energy analysis interanalyzer.py
# Usage: python Interface Analyzer 
########################################################################
# import all the necessary functions required to make functional calls
import  sys, re, os

amino_acids = ['A', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
#amino_acids = ['A', 'T', 'Y']

if len(sys.argv)<5:
     print('requires at least 4 arguments %s given' % (len(sys.argv)))
     print('Usage for ddg resfile format: make_mut_file.py <number of mutations> <aa_input_pdb> <res_number> <aa_2_mutate_2>')
     print('Usage for canonical resfile format : make_mut_file.py <aa_scan> <aa_input_pdb> <res_number :{PDBNUM[<ICODE>]|POSENUM}> <chain id> <command:PIKAA> ')
     print('PDBNUM; pdb numbering 5 or with icode 5A; POSENUM; rosetta renumbered pdb or from another rosetta application')
     print('command : controls sequence identity see rosetta control acronyms PIKAA requires specification of residue to change to either one or several. all other commands do not require specification so aa_2_mutate_2 should be left out') 
     sys.exit(0)
if sys.argv[4]!='PIKAA':
  n=sys.argv[1]
  wt=sys.argv[2]
  res_number=sys.argv[3]
  mut=sys.argv[4]
  p=open('%s%s%s.file'%(wt,res_number,mut),'w')
  p.write('total %s\n%s\n%s %s %s' %(n,n,wt,res_number,mut))
  p.close()
if sys.argv[4]=='PIKAA':
  wt=sys.argv[1]
  res_number=sys.argv[2]
  chain=sys.argv[3]
  command=sys.argv[4]
  for res in amino_acids:
     if res==wt:continue
     p=open('%s%s%s.resfile'%(wt,res_number,res),'w')
     p.write('NATAA\nSTART\n\n%s %s %s %s' %(res_number,chain,command,res))
     p.close()




