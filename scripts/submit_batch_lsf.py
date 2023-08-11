#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
---------------------------------------------------------------
 authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
          H. Garcia-Tecocoatzi
---------------------------------------------------------------
"""
# lsf job submitter
# usage: python3 scripts/submit_batch_lsf.py three_quark omegas

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

# prepare lsf option file
classad='''
#BSUB -J quark[1-{4}]
#BSUB -L /bin/bash
#BSUB -n 1
#BSUB -o {1}/output_batch/{3}/%I.out
#BSUB -e {1}/output_batch/{3}/%I.err
{0} {1}/{2} $LSB_JOBINDEX {1} {3}
'''.format(py3_path, workpath, run_bottom, baryons, n_jobs)

logpath = '.'

with open(logpath+'/lsf_options.jdl','w') as jdlfile:
  jdlfile.write(classad)

print("************************ Batch jobs for: ", baryons, " masses and decays  ************************")
system('bsub -q long < %s/lsf_options.jdl'%logpath)
