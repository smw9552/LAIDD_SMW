# -*- coding: utf-8 -*-
# Author: Myungwon Seo
# Date: 2023-10-20
# E-mail: mwseo@krict.re.kr; seomyungwon@gmail.com

import pandas as pd
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys

#FP list: MACCS, RDKit, Morgan/Circular

# 파일경로 설정
current_directory = os.path.dirname(__file__)
input_file = os.path.join(current_directory, "..", "dataset", "NPASS_NPs.xlsx")

# 파일경로 직접 설정
#input_file = "userpath\\NPASS_NPs.xlsx"
input_file = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\NPASS_NPs.xlsx"
output_file1 = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\output\\RDKit_output.xlsx"
output_file2 = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\output\\RDKit_output.xlsx"
output_file3 = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\output\\RDKit_output.xlsx"

#파일을 dataframe 형태로 불러오기
df = pd.read_excel(input_file)

#파일속성 확인
df.columns

#SMILES 정보 추출하기
smiles_list = df['Smiles']
cas_list = df['CAS']
cid_list = df['PubChem_CID']
mol_list = []

MACCS_list = []
RDKit_list = []
Morgan_list = []


#Mol file 생성
for ai in range(0, len(smiles_list)):
    mol_list.append(Chem.MolFromSmiles(smiles_list[ai]))


#generate bit column (파일 작성을 위해 bit column 생성)
temp_smiles = 'CCO'
temp_mol = Chem.MolFromSmiles(temp_smiles)
temp_MACCS = MACCSkeys.GenMACCSKeys(temp_mol) #167
temp_RDKit = Chem.RDKFingerprint(temp_mol) # 2048
temp_Morgan = AllChem.GetMorganFingerprintAsBitVect(temp_mol, radius = 6, nBits=1024) #1024 (default)
temp_MACCS_binary = temp_MACCS.ToBitString()
temp_RDKit_binary = temp_RDKit.ToBitString()
temp_Morgan_binary = temp_Morgan.ToBitString()

for a in range(0, len(temp_MACCS_binary)):
    column = "MACCS" + "_" + str(a+1)
    MACCS_list.append(column)

for b in range(0, len(temp_RDKit_binary)):
    column = "RDKit" + "_" + str(b+1)
    RDKit_list.append(column)

for c in range(0, len(temp_Morgan_binary)):
    column = "ECFP6" + "_" + str(c + 1)
    Morgan_list.append(column)


for bi in range(0, len(mol_list)):
    print("count:", bi)
    print(smiles_list[bi])

    # 구조 생성문제 확인하기 위한 작업
    type_check = str(mol_list[bi])

    if (type_check == 'None'):
        # Smile 정보로 mol 생성이 불가능한 경우 descriptor 계산이 안되므로 모두 0으로 처리 (에러방지)
        # 0 정보를 담을 임시 리스트 작성
        temp_MACCS_FP = []
        temp_RDKit_FP = []
        temp_Morgan_FP = []

        for bii_1 in range(0, 167):
            temp_MACCS_FP.append(0)
        for bii_2 in range(0, 2048):
            temp_RDKit_FP.append(0)
        for bii_3 in range(0, 1024):
            temp_Morgan_FP.append(0)

        # 0 정보가 담긴 리스트를 추가
        MACCS_FP.append(temp_MACCS_FP)
        RDKit_FP.append(temp_RDKit_FP)
        Morgan_FP.append(temp_Morgan_FP)

    else:
        MACCS = MACCSkeys.GenMACCSKeys(mol_list[bi])  # 167
        RDKit = Chem.RDKFingerprint(mol_list[bi])  # 2048
        Morgan = AllChem.GetMorganFingerprintAsBitVect(mol_list[bi], radius=6, nBits=1024)  # 1024 (default)
        MACCS_binary = MACCS.ToBitString()
        RDKit_binary = RDKit.ToBitString()
        Morgan_binary = Morgan.ToBitString()

        # 각 분자별 FP 최종정리 (리스트 내부의 리스트를 풀어서 하나의 리스트로 작성)
        final_MACCS_FP = []
        final_RDKit_FP = []
        final_Morgan_FP = []

        for item_maccs in MACCS_binary:
            if isinstance(item_maccs, list):
                final_MACCS_FP.extend(item_maccs)
            else:
                final_MACCS_FP.append(item_maccs)
        MACCS_FP.append(final_MACCS_FP)

        for item_rdkit in RDKit_binary:
            if isinstance(item_rdkit, list):
                final_RDKit_FP.extend(item_rdkit)
            else:
                final_RDKit_FP.append(item_rdkit)
        RDKit_FP.append(final_RDKit_FP)

        for item_morgan in Morgan_binary:
            if isinstance(item_morgan, list):
                final_Morgan_FP.extend(item_morgan)
            else:
                final_Morgan_FP.append(item_morgan)
        Morgan_FP.append(final_Morgan_FP)

result_MACCS = np.array(MACCS_FP).astype(int)
result_RDKit = np.array(RDKit_FP).astype(int)
result_Morgan = np.array(Morgan_FP).astype(int)

MACCS_df = pd.DataFrame(data=result_MACCS, columns=MACCS_list)
RDKit_df = pd.DataFrame(data=result_RDKit, columns=RDKit_list)
Morgan_df = pd.DataFrame(data=result_Morgan, columns=Morgan_list)

# df에 CID, label column 추가
# MACCS_df.insert(loc=0, column="CID", value=cid_list)
MACCS_df.insert(loc=0, column="CAS", value=cas_list)
MACCS_df.insert(loc=1, column="PubChem_CID", value=cid_list)
MACCS_df.to_excel(output_file_1, index=False)

# RDKit_df.insert(loc=0, column="CID", value=cid_list)
RDKit_df.insert(loc=0, column="CAS", value=cas_list)
RDKit_df.insert(loc=1, column="PubChem_CID", value=cid_list)
RDKit_df.to_excel(output_file_2, index=False)

# Morgan_df.insert(loc=0, column="CID", value=cid_list)
Morgan_df.insert(loc=0, column="CAS", value=cas_list)
Morgan_df.insert(loc=1, column="PubChem_CID", value=cid_list)
Morgan_df.to_excel(output_file_3, index=False)


