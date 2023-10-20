from rdkit import Chem
import pandas as pd
import os
from mordred import Calculator, descriptors


# 파일경로 설정
current_directory = os.path.dirname(__file__)
input_file = os.path.join(current_directory, "..", "dataset", "NPASS_NPs.xlsx")


# 파일경로 직접 설정
#input_file = "userpath\\NPASS_NPs.xlsx"
input_file = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\NPASS_NPs.xlsx"
output_file = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\output\\Mordred_output.xlsx"

#파일을 dataframe 형태로 불러오기
df = pd.read_excel(input_file)

#파일속성 확인
df.columns

#SMILES 정보 추출하기
smiles_list = df['Smiles']
cas_list = df['CAS']
cid_list = df['PubChem_CID']
mol_list = []


#Mol file 생성
for ai in range(0, len(smiles_list)):
    mol_list.append(Chem.MolFromSmiles(smiles_list[ai]))



remove_processing = False
replace_precessing = True # 1,613개의 descriptor 수를 고정하기 위한 방법

#Descriptor calculator 생성
calc = Calculator(descriptors, ignore_3D=True)

# 일부만 선택하여 구동
# selected_descriptor = [descriptors.ABCIndex]

# calculate descriptors as pandas (descriptor 전체 구동)
df_desc = calc.pandas(mol_list)

if (remove_processing == True):
    # chage data type and preprocessing
    del_index = []
    for bi in range(0, len(df_desc.columns)):
        if (df_desc.dtypes[bi] == 'object' or df_desc.dtypes[bi] == 'bool'):
            del_index.append(bi)

    df_desc = df_desc.drop(df_desc.columns[del_index], axis=1)  # 문자가 포함된 열 제거
    # df_desc = df_desc.astype('int64') #type 고정

if (replace_precessing == True):  # 계산되는 column 숫자를 1,613개로 통일시키기 위해 계산이 안되거나 문자가 출력되는 부분을 0으로 치환
    rep_index = []
    for ci in range(0, len(df_desc.columns)):
        if (df_desc.dtypes[ci] == 'object' or df_desc.dtypes[ci] == 'bool'):
            rep_index.append(ci)

    for di in range(0, len(rep_index)):
        df_desc.iloc[:, rep_index[di]] = 0

df_desc.insert(0, "CAS", cas_list)
df_desc.insert(1, "PubChem_CID", cid_list)

# write output file
df_desc.to_excel(output_file, index=False)