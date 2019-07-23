import uproot
import pandas as pd
import numpy as np
import awkward

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_nano(infile):
  uf = uproot.open(infile)
  tt = uf['Events']
  arrs = tt.arrays(tt.keys())
  ## vvs = ['pt', 'eta', 'phi']
  ## els = awkward.JaggedArray.zip(  *[arrs['Electron_%s' % i] for i in vvs])
  ## muons = awkward.JaggedArray.zip(*[arrs['Muon_%s'  % i] for i in vvs])
  ## #trks = awkward.JaggedArray.zip( *[arrs['Track_%s' % i] for i in vvs])
  ## arrs['muons'] = muons
  ## arrs['eles'] = els
  ## arrs['BToKEE'] = awkward.JaggedArray.zip(
  ##   arrs['BToKEE_pt'],
  ##   arrs['BToKEE_eta'],
  ##   arrs['BToKEE_phi'],
  ##   arrs['BToKEE_mass'],
  ##   arrs['BToKEE_svprob'],
  ##   arrs['BToKEE_cos2D'],
  ##   arrs['BToKEE_mll'],
  ##   )
  ## arrs['BToKMuMu'] = awkward.JaggedArray.zip(
  ##   arrs['BToKMM_pt'],
  ##   arrs['BToKMM_eta'],
  ##   arrs['BToKMM_phi'],
  ##   arrs['BToKMM_mass'],
  ##   arrs['BToKMM_svprob'],
  ##   arrs['BToKMM_cos2D'],
  ##   arrs['BToKMM_mll'],
  ##   )
  ## #arrs['trks'] = trks
  return arrs


def byval_validation(v1, v2):
  return (v1 == v2).all()

def stat_validation(v1, v2, name = '', nbins = 20):
  if v1.shape[0] == 0 and v2.shape[0] == 0:
    return True
  elif v1.shape[0] == 0 or v2.shape[0] == 0:
    return False
  M = max(v1.max(), v2.max())*1.2
  m = min(v1.min(), v2.min())*0.9
  plt.clf()
  h1, _, _ = plt.hist(v1, range = (m,M), bins = nbins, label = 'old', histtype = 'step')
  h2, _, _ = plt.hist(v2, range = (m,M), bins = nbins, label = 'new', histtype = 'step')
  plt.legend(loc='best')
  plt.savefig('validation/%s.png' % name)
  return (h1 == h2).all()

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('f_old', help='file path')
parser.add_argument('f_new', help='file path')
args = parser.parse_args()

old = load_nano(args.f_old)
new = load_nano(args.f_new)

old_k = set(old.keys())
new_k = set(new.keys())
intersection = old_k.intersection(new_k)

if len(intersection) != len(old_k) or len(intersection) != len(new_k):
  raw_input("The branches are different!")
  
for branch in intersection:
  v_old = old[branch]
  v_new = new[branch]
  if hasattr(v_old, 'flatten'):
    v_old = v_old.flatten()
    v_new = v_new.flatten()
  valid = stat_validation(v_old, v_new, branch)
  if valid:
    print branch, '--> OK!'
  else:
    print branch, '--> FAIL!'

matching = {
  'BToKEE_kIdx'   : 'Track_',
  'BToKEE_el1Idx' : 'Electron_',
  'BToKEE_el2Idx' : 'Electron_',
  
  'BToKMuMu_kIdx'   : 'Track_',
  'BToKMuMu_mu1Idx' : 'Muon_',
  'BToKMuMu_mu2Idx' : 'Muon_',
}

for key, match in matching.iteritems():
  for branch in intersection:
    if not branch.startswith(match): continue
    v_old = old[branch][old[key]].flatten()
    v_new = new[branch][new[key]].flatten()
    name = '_'.join([key, branch])
    valid = stat_validation(v_old, v_new, name)
    if valid:
      print name, '--> OK!'
    else:
      print name, '--> FAIL!'
