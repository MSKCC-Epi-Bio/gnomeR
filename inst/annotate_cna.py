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

def cnaAnnotate(in_cna, out_cna, clin_file, token):
  setoncokbapitoken(token)
  getcuratedgenes()

  default_cancer_type = 'cancer'
  cancertypemap = {}
  if clin_file != '':
    readCancerTypes(clin_file, cancertypemap)
    
  log.info('annotating CNAs %s ...' % in_cna)
  processcnagisticdata(in_cna, out_cna, False, default_cancer_type,
                         cancertypemap, True, False)
  # argv.input_file, argv.output_file, argv.previous_result_file, argv.default_cancer_type,
  #                        cancertypemap, True, argv.annotate_gain_loss
  log.info('CNAs done!')
