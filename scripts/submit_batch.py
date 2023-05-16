#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
---------------------------------------------------------------
 authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          H. Garcia-Tecocoatzi
---------------------------------------------------------------
"""
# htcondor job submitter
# usage: python3 scripts/submit_batch.py three_quark/diquark

import sys
import json
import os
from os import system,getenv,getuid,getcwd,popen

workpath = getcwd()

if len(sys.argv)!=3:
  sys.exit('Please indicate wich quark structure to run for batch jobs')

baryons = "omegas"
n_jobs  = 10
with open(workpath+"/config/three_quark_config.json", "r") as jsonfile:
    config = json.load(jsonfile)
# baryons = config["baryons"] bug happening, need to feed directly the baryon name
n_jobs  = config["n_jobs"]

three_di = sys.argv[1]
baryons = sys.argv[2]
py3_path = popen('which python3').read().strip()

if (three_di == 'three_quark'):
  run_bottom = 'scripts/bootstrap_three_quark.py'
elif (three_di=='diquark'):
  run_bottom = 'bottom_bootstrap_diquark_flavor.py'
else:
  sys.exit('quark structure not supported')

classad='''
universe = vanilla
executable = {0}
getenv = True
arguments = {1}/{2} $(Process) {1} {3}
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
