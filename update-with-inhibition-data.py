#!/usr/bin/python3
import os
import sys
import csv
import pandas
import requests
import subprocess
from rdkit import Chem

def inchi_string(smiles):
    smiles = smiles.replace('+H', 'H+') # annoying MOE variant RDKit doesn't like
    return Chem.MolToInchi(Chem.MolFromSmiles(smiles))

def inchi_isomer(smiles):
    return '/'.join(inchi_string(smiles).split('/')[:4])

IDX = -1
def primary_key(row):
    global IDX
    IDX += 1
    return IDX

def best_match(smiles):
    inchi = inchi_string(smiles)
    found = invitro_data['INCHI'] == inchi
    indices = list(found[found == True].index)
    if not indices:
        # fallback to fuzzy
        fuzzy = inchi_isomer(smiles)
        found = invitro_data['FUZZY'] == fuzzy
        indices = list(found[found == True].index)

    if not indices:
        return -1

    if len(indices) > 1:
        print("Multiple matches found for", smiles)

    return indices[0]

invitro_data = pandas.read_csv('doi.org-10.1101-2020.04.03.023846-with-smiles.csv')
invitro_data['INCHI'] = invitro_data['SMILES'].apply(inchi_string)
invitro_data['FUZZY'] = invitro_data['SMILES'].apply(inchi_isomer)
invitro_data['IDX'] = invitro_data.apply(primary_key, axis=1)

data_to_merge = pandas.read_csv(sys.argv[1])
data_to_merge['MERGE_IDX'] = data_to_merge['SMILES'].apply(best_match)

merged = data_to_merge.merge(invitro_data, how='inner', left_on=['MERGE_IDX'], right_on=['IDX'])
merged.to_csv(os.path.splitext(os.path.basename(sys.argv[1]))[0] + '-with-inhibition-index.csv')

print('Found', len(merged), 'matches')

from matplotlib import pyplot as plt
plt.xlabel('Inhibition Index')
plt.ylabel(sys.argv[2])
plt.scatter(merged['Inhibition Index'], merged[sys.argv[2]])
plt.savefig('Inhibition Index vs ' + sys.argv[2] + '.png')

