from rdkit.Chem import Descriptors
from rdkit import Chem
import pandas as pd
import numpy as np
import os
from padelpy import from_smiles


# 파일경로 설정
current_directory = os.path.dirname(__file__)
input_file = os.path.join(current_directory, "..", "dataset", "NPASS_NPs.xlsx")


# 파일경로 직접 설정
#input_file = "userpath\\NPASS_NPs.xlsx"
input_file = "C:\\Users\\user\\PycharmProjects\\LAIDD_SMW\\dataset\\NPASS_NPs.xlsx"
output_file_1 = "C:\\Users\\user\\PycharmProjects\\LAIDD_SMW\\dataset\\output\\PaDelPy_desc_output.xlsx"
output_file_2 = "C:\\Users\\user\\PycharmProjects\\LAIDD_SMW\\dataset\\output\\PaDelPy_fp_output.xlsx"

#파일을 dataframe 형태로 불러오기
df = pd.read_excel(input_file)

#파일속성 확인
df.columns

#SMILES 정보 추출하기
smiles_list = df['Smiles']
smiles_list = list(smiles_list) # descriptor 계산을 위한 리스트 형태로 변환
cas_list = df['CAS']
cid_list = df['PubChem_CID']
mol_list = []

# calculate molecular descriptors
desc = from_smiles(smiles_list, fingerprints=False, descriptors=True)


# (descriptor) dataframe으로 type 변환
temp_df = pd.DataFrame()
merge_df = pd.DataFrame()
for bi in range(0, len(desc)):

    temp_df = pd.DataFrame([desc[bi]])
    merge_df = pd.concat([temp_df, merge_df])
    print(merge_df)

#빈값을 0으로 전환
merge_df = merge_df.replace('',0)


# calculate PubChem fingerprint
desc_fp = from_smiles(smiles_list, fingerprints=True, descriptors=False)
desc_fp = pd.Series(desc_fp)

# (molecular fingerprint) dataframe으로 type 변환
temp_fp_df = pd.DataFrame()
merge_fp_df = pd.DataFrame()
for ci in range(0, len(desc_fp)):

    temp_fp_df = pd.DataFrame([desc_fp[ci]])
    merge_fp_df = pd.concat([temp_fp_df, merge_fp_df])
    print(merge_fp_df)

#빈값을 0으로 전환
merge_fp_df = merge_fp_df.replace('',0)

# insert info as new column
merge_df.insert(0, "CAS", cas_list)
merge_df.insert(1, "PubChem_CID", cid_list)

merge_fp_df.insert(0, "CAS", cas_list)
merge_fp_df.insert(1, "PubChem_CID", cid_list)



#save descriptors
merge_df.to_excel(output_file_1, index = False)
merge_fp_df.to_excel(output_file_1, index = False)