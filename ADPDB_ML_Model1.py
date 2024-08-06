#Sai
#Sai
#OM Sai Ram

import pickle
from peptides import Peptide as pept
from Bio.SeqUtils.ProtParam import ProteinAnalysis as proana
import streamlit as lit
from itertools import product
import pandas as pd
## Models
RandomForestClassifier = pickle.load(open('RandomForestClassifier.pkl', 'rb'))
RF_scaler = pickle.load(open('RandomForestClassifier_scaler.pkl', 'rb'))
ExtraTreesClassifier = pickle.load(open('ExtraTreesClassifier.pkl', 'rb'))
ET_scaler = pickle.load(open('ExtraTreesClassifier_scaler.pkl', 'rb'))
## Functions
def classifier(name,inputs):
    match name:
        case 'RandomForestClassifier':
            _ = [inputs.pop(i) for i in [20, 21, 22, 23, 24, 25]]#
            inputs = RF_scaler.transform([inputs])
            out = RandomForestClassifier.predict(inputs)
            prob = RandomForestClassifier.predict_proba(inputs)
        case 'ExtraTreesClassifier':
            _ = [inputs.pop(i) for i in [21,23,24]]#
            inputs = ET_scaler.transform([inputs])
            out = ExtraTreesClassifier.predict(inputs)
            prob = ExtraTreesClassifier.predict_proba(inputs)
    return out.max(),prob.max()

def descriptor(seq):
    pep = pept(seq.strip())
    proan = proana(seq.strip())
    _ = (
            list(pep.frequencies().values()) + 
            [pep.aliphatic_index()] + 
            [pep.instability_index()] + 
            [pep.hydrophobicity()] + 
            [pep.hydrophobic_moment()] + 
            [pep.isoelectric_point()] + 
            [pep.molecular_weight()] + 
            [pep.charge()] +
            [proan.aromaticity()]
        )
    return _
## UI
lit.set_page_config(layout='wide')
lit.write("""
# Welcome to the ADPDB Sequence Prediction Tool!
*A tool for prediction of anti-dengue peptides.*
""")
sequence = lit.text_area('Please enter the input sequence (Plain text/FASTA/multi-FASTA formats supported)')
file_query = lit.file_uploader("Or, you may upload file (Plain text/FASTA/multi-FASTA formats supported)")
if file_query:
    sequence = StringIO(file_query.getvalue().decode("utf-8")).read().upper()
repeat = 0
#if lit.checkbox('Do you want list of permutation of all input amino acids?'):
#    repeat = lit.text_input('What length of peptides should be created? (Please enter integer value)')
model = lit.radio('Choose a model:', ['RandomForestClassifier', 'ExtraTreesClassifier'])
submit = lit.button('Submit')
if submit:
    if sequence:
        if sequence.count('>') == 1:
            sequence = ''.join(sequence.split('\n')[1:])
        elif sequence.count('>') > 1:
            sequence = [''.join(i.split('\n')[1:]) for i in sequence.split('>') if ' ' not in i]
        if repeat:
            repeat = int(repeat)
            sequence = [''.join(i) for i in product(sequence, repeat=repeat)]
        if isinstance(sequence, str):
            sequence = [sequence]
        df = []
        for i in sequence:
            descriptors = descriptor(i)
            out,prob = classifier(model, descriptors)
            df.append([i, ('Anti-Dengue' if out==1 else'Non Anti-Dengue'), prob])
        df = pd.DataFrame(df, columns=['Sequence', 'Class', 'Probability '])
        df = df.sort_values(by=['Probability'], ascending=False)
        df.reset_index(drop=True, inplace=True)
        lit.dataframe(df)
    else:
        lit.warning('NO input given!')
