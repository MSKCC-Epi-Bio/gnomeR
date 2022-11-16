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

def mafAnnotate(in_maf, out_maf, clin_file, token):
  # if argv.sample_ids_filter:
  #   setsampleidsfileterfile(argv.sample_ids_filter)
  # if argv.cancer_hotspots_base_url:
  #   setcancerhotspotsbaseurl(argv.cancer_hotspots_base_url)
  # if argv.oncokb_api_url:
  # oncokbapiurl = "https://www.oncokb.org/api/v1"
  # setoncokbbaseurl(oncokbapiurl)
  setoncokbapitoken(token)
  # getcuratedgenes()

  default_reference_genome = None
  default_cancer_type = 'cancer'
  cancertypemap = {}
  if clin_file != '':
    readCancerTypes(clin_file, cancertypemap)
    # default_cancer_type = False
    
  # log.info('CancerMap %s ...' % cancertypemap)
  # log.info('default_cancer_type %s ...' % default_cancer_type)
  
  log.info('annotating MAF %s ...' % in_maf)
  processalterationevents(in_maf, out_maf, False, default_cancer_type, cancertypemap, True, False, default_reference_genome)
#argv.input_file, argv.output_file, argv.previous_result_file, argv.default_cancer_type, cancertypemap, True, argv.annotate_hotspot
  log.info('MAF done!')
