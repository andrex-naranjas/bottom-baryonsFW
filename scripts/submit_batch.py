#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
---------------------------------------------------------------
 authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          H. Garcia-Tecocoatzi
 ---------------------------------------------------------------
'''

# htcondor job submitter
# python3 scripts/submit_batch.py baryons three_quark/diquark n_jobs

import sys
from os import system,getenv,getuid,getcwd,popen

workpath = getcwd()

if len(sys.argv)!=4:
  sys.exit('Not enough arguments to set the batch job')

baryons  = sys.argv[1]
three_di = sys.argv[2]
n_jobs = sys.argv[3]
py3_path = popen('which python3').read().strip()

if three_di=='three_quark':
  run_bottom = 'scripts/bootstrap_three_quark.py'
elif three_di=='diquark':
  run_bottom = 'bottom_bootstrap_diquark.py'
else:
  sys.exit('quark structure not supported')

classad='''
universe = vanilla
executable = {0}
getenv = True
arguments = {1}/{2} {3} $(Process) {1}
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = {1}/output_batch/{3}/$(Process).out
error = {1}/output_batch/{3}/$(Process).err
log = {1}/output_batch/{3}/$(Process).log
Queue {4}
'''.format(py3_path, workpath, run_bottom, baryons, n_jobs)

logpath = '.'

with open(logpath+'/condor.jdl','w') as jdlfile:
  jdlfile.write(classad)

print("************************ Batch jobs for: ", baryons, " masses and decays  ************************")
#print('condor_submit %s/condor.jdl'%logpath)
system('condor_submit %s/condor.jdl'%logpath)
