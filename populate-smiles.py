#!/usr/bin/python3

import csv
import pandas
import requests
import subprocess
from rdkit import Chem

dtable = pandas.read_csv('doi.org-10.1101-2020.04.03.023846.csv')

#import psycopg2 as pg
#conn = pg.connect(... postgres connection string ...)

import sqlite3
conn = sqlite3.connect('chembl_26/chembl_26_sqlite/chembl_26.db')

cursor = conn.cursor()

SQL = """ 
select 
compound_structures.canonical_smiles as compound_smiles, 
molecule_dictionary.chembl_id as compound_chembl_id, 
molecule_dictionary.pref_name as compound_name, 
molecule_dictionary.max_phase as drug_development_phase 
FROM compound_structures 
JOIN molecule_dictionary 
ON compound_structures.molregno = molecule_dictionary.molregno """

def search_chembl_for_compound(chemical_name):
    if chemical_name.startswith('CHEMBL'):
        sql = SQL + "WHERE molecule_dictionary.chembl_id=?;"
    else:
        sql = SQL + "WHERE molecule_dictionary.pref_name LIKE ? COLLATE NOCASE;"

    cursor.execute(sql, (chemical_name,))
    records = list(cursor)
    if records:
        return records

    for part in chemical_name.split():
        cursor.execute(sql, (part,))
        records = list(cursor)
        if records:
            return records

    for part in chemical_name.split('-'):
        cursor.execute(sql, (part,))
        records = list(cursor)
        if records:
            return records

    return []


UNKNOWN_TO_CHEMBL = csv.writer(open('unknown-to-chembl.csv', 'w'))

def query_nih_resolver(chemical_name):
    for query_str in ([chemical_name] + chemical_name.split()):
        response = requests.get(f'http://cactus.nci.nih.gov/chemical/structure/{query_str}/stdinchi')
        if response.status_code == 200: 
            return query_str, response.text

    print(f"NIH can not resolve {chemical_name}")
    return None, None

def search_nih(chemical_name):
    query_str, inchi = query_nih_resolver(chemical_name)
    if query_str is None:
        return []

    sql = SQL + "WHERE compound_structures.standard_inchi LIKE ?;"
    cursor.execute(sql, (inchi,))
    records = list(cursor)
    if not records:
        response = requests.get(f'http://cactus.nci.nih.gov/chemical/structure/{query_str}/smiles')
        print(response.text, query_str, "is not in ChEMBL")
        UNKNOWN_TO_CHEMBL.writerow((response.text, query_str, chemical_name))

    return records

def inchi_isomer(smiles):
    return '/'.join(Chem.MolToInchi(Chem.MolFromSmiles(smiles)).split('/')[:4])

# constructed by manual google search
RENAME_TABLE = {
    'Nitrofural' : 'NITROFURAZONE',
    'Zoledronic acid hydrate' : 'ZOLEDRONIC ACID',
    'Dehydroisoandosterone 3-acetate' : 'PRASTERONE',
    'Fluspirilen' : 'FLUSPIRILENE',
    'Cresopirine' : 'CHEMBL1234172',
    'Rifampicin' : 'RIFAMPIN',
    'Azacytidine-5' : 'AZACITIDINE',
    'Eserine hemisulfate salt' : 'PHYSOSTIGMINE',
    'Scopolamin-N-oxide hydrobromide' : 'Too many isomers in ChEMBL to choose from',
    'Acetylsalicylsalicylic' : 'CHEMBL350343',
    'Oxibendazol' : 'CHEMBL1087630',
    'Gentamicine' : 'CHEMBL463809',
}

def find_compound(chemical_name):
    if chemical_name in RENAME_TABLE:
        chemical_name = RENAME_TABLE[chemical_name]
    records = search_chembl_for_compound(chemical_name)
    if not records:
        records = search_nih(chemical_name)

    if not records:
        return None
    elif len(records) > 1:
        all_smiles = set()
        all_chemblid = set()
        all_compound_name = set()
        all_max_phase = set()
        
        for smiles, chemblid, compound_name, max_phase in records:
            longest_smi = list(sorted(smiles.split('.'), key=len, reverse=True))[0]
            all_smiles.add(longest_smi)
            all_chemblid.add(chemblid)
            all_compound_name.add(compound_name)
            all_max_phase.add(max_phase)

        if len(all_smiles) > 1 and len(set(inchi_isomer(s) for s in all_smiles)) > 1:
            print(f"Multiple records found for {chemical_name}:", records)

        records = [(all_smiles.pop(), ','.join(sorted(all_chemblid)), ','.join(sorted(all_compound_name)), ','.join(map(str, sorted(all_max_phase))))]

    return records[0]


with open('doi.org-10.1101-2020.04.03.023846-with-smiles.csv', 'w') as ofs:
    writer = csv.writer(ofs)
    writer.writerow(('SMILES', 'CHEMBL_ID', 'CHEMBL_COMPOUND_NAME', 'CHEMBL_MAX_PHASE',) + tuple(map(str, dtable.columns)))

    nfound = 0
    for chemical_name in dtable['Chemical name']:
        compound = find_compound(chemical_name)
        if compound is None:
            print(f"No record found for {chemical_name}")
            continue

        nfound += 1
        values = dtable[dtable['Chemical name'] == chemical_name].values
        if len(values) > 1:
            print(f"Multiple inhibitions found for {chemical_name}:", values)
        orig_data = tuple(map(str, values[0]))
        writer.writerow(compound + orig_data)

    print(nfound, 'SMILES found for',  len(dtable), 'data points')
