import argparse
# from AnnotatorCore import *
import sys
import csv
import requests
import os.path
import logging
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import date
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger('OncoKB Annotator')

def fusionAnnotate(in_fusion, out_fusion, clin_file, token):
  setoncokbapitoken(token)
  # getcuratedgenes()

  default_cancer_type = 'cancer'
  cancertypemap = {}
  if clin_file != '':
    readCancerTypes(clin_file, cancertypemap)
    
  log.info('annotating fusions %s ...' % in_fusion)
  processsv(in_fusion, out_fusion, False, default_cancer_type, cancertypemap, None)
  # processsv(argv.input_file, argv.output_file, argv.previous_result_file, argv.default_cancer_type,
  #             cancertypemap, argv.structural_variant_name_format)
  log.info('Fusions done!')
  
  # processalterationevents(in_maf, out_maf, False, default_cancer_type, cancertypemap, True, False)
