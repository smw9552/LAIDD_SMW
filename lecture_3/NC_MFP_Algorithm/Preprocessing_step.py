# -*- coding: utf-8 -*-
# Author: Myungwon Seo
# Date: 2023-10-20
# E-mail: mwseo@krict.re.kr; seomyungwon@gmail.com


from rdkit import Chem
import pandas as pd

class Preprocessing:


    def get_Query_Mol_single(self):

        #Query molecule name: Vismiaguianin A
        qMol = Chem.MolFromSmarts('[#8]1-[#6](=[#6](-[#6](=[#8])-[#8]2)-[#6](-[#1])=[#6](-[#1])-[#6]1(-[#6](-[#1])(-[#1])-[#1])-[#6](-[#1])(-[#1])-[#1])-[c]3:[c]2:[c]4:[c](:[c](-[#1]):[c](:[c](-[#1]):[c]4-[#8]-[#1])-[#8]-[#6](-[#1])(-[#1])-[#1]):[c](-[#1]):[c]3-[#6](-[#1])(-[#1])-[#1]')

        return qMol


    def get_Query_Mols(self, FilePath):

        f = open(FilePath, 'r')
        qMols = []
        while True:
            line = f.readline()
            if not line:break
            new_line = line.split("\t")
            qMols.append(Chem.MolFromSmarts(new_line[1].strip()))

        return qMols


    def get_Query_Mol_single_Smarts(self):

        #Query molecule name: Vismiaguianin A
        qMol_Smarts = '[#8]1-[#6](=[#6](-[#6](=[#8])-[#8]2)-[#6](-[#1])=[#6](-[#1])-[#6]1(-[#6](-[#1])(-[#1])-[#1])-[#6](-[#1])(-[#1])-[#1])-[c]3:[c]2:[c]4:[c](:[c](-[#1]):[c](:[c](-[#1]):[c]4-[#8]-[#1])-[#8]-[#6](-[#1])(-[#1])-[#1]):[c](-[#1]):[c]3-[#6](-[#1])(-[#1])-[#1]'

        return qMol_Smarts

    def get_Query_Mols_Smarts(self, FilePath):

        f = open(FilePath, 'r')
        qMols_Smarts = []
        while True:
            line = f.readline()
            if not line:break
            new_line = line.split("\t")
            qMols_Smarts.append(new_line[1].strip())

        return qMols_Smarts


    def get_Query_Mols_Smarts_new(self, FilePath):

        #smiles, cas, cid 정보가 있는 경우 활용 가능 (smiles -> mol -> smarts)
        df = pd.read_excel(FilePath)
        smiles_list = df['Smiles']
        cas_list = df['CAS']
        cid_list = df['PubChem_CID']
        mol_list = []
        smart_list = []

        # Mol file 생성
        for ai in range(0, len(smiles_list)):
            mol = Chem.MolFromSmiles(smiles_list[ai])
            mol_with_H = Chem.AddHs(mol)
            mol_list.append(mol_with_H)

        for bi in range(0, len(mol_list)):
            smart_list.append(Chem.MolToSmarts(mol_list[bi]))

        return smart_list







