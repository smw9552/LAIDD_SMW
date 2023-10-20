# -*- coding: utf-8 -*-
# Author: Myungwon Seo
# Date: 2023-10-20
# E-mail: mwseo@krict.re.kr; seomyungwon@gmail.com

from rdkit.Chem import Descriptors
from rdkit import Chem
import pandas as pd
import numpy as np
import os


# 파일경로 설정
current_directory = os.path.dirname(__file__)
input_file = os.path.join(current_directory, "..", "dataset", "NPASS_NPs.xlsx")


#Single_descriptor name (106개)
rdkit_desc_list = ["BalabanJ", "BertzCT", "Ipc", "HallKierAlpha", "Kappa1", "Kappa2", "Kappa3", "Phi", "Chi0", "Chi1",
                     "Chi0n", "Chi1n", "Chi2n", "Chi3n", "Chi4n", "Chi0v", "Chi1v", "Chi2v", "Chi3v", "Chi4v", "MolLogP",
                     "MolMR", "MolWt", "ExactMolWt", "HeavyAtomCount", "HeavyAtomMolWt", "NHOHCount", "NOCount", "NumHAcceptors",
                     "NumHDonors", "NumHeteroatoms", "NumRotatableBonds", "NumValenceElectrons", "NumAmideBonds", "NumAromaticRings",
                     "NumSaturatedRings", "NumAliphaticRings", "NumAromaticCarbocycles", "NumAromaticHeterocycles", "NumSaturatedCarbocycles",
                     "NumSaturatedHeterocycles", "NumAliphaticCarbocycles", "NumAliphaticHeterocycles", "RingCount", "FractionCSP3",
                     "NumSpiroAtoms", "NumBridgeheadAtoms", "TPSA", "LabuteASA", "PEOE_VSA1", "PEOE_VSA2", "PEOE_VSA3", "PEOE_VSA4", "PEOE_VSA5",
                     "PEOE_VSA6", "PEOE_VSA7", "PEOE_VSA8", "PEOE_VSA9", "PEOE_VSA10", "PEOE_VSA11", "PEOE_VSA12", "PEOE_VSA13", "PEOE_VSA14",
                     "SMR_VSA1", "SMR_VSA2", "SMR_VSA3", "SMR_VSA4", "SMR_VSA5", "SMR_VSA6", "SMR_VSA7", "SMR_VSA8", "SMR_VSA9", "SMR_VSA10",
                     "SlogP_VSA1", "SlogP_VSA2", "SlogP_VSA3", "SlogP_VSA4", "SlogP_VSA5", "SlogP_VSA6", "SlogP_VSA7", "SlogP_VSA8", "SlogP_VSA9",
                     "SlogP_VSA10", "SlogP_VSA11", "SlogP_VSA12", "EState_VSA1", "EState_VSA2", "EState_VSA3", "EState_VSA4", "EState_VSA5",
                     "EState_VSA6", "EState_VSA7", "EState_VSA8", "EState_VSA9", "EState_VSA10", "EState_VSA11", "VSA_EState1", "VSA_EState2",
                     "VSA_EState3", "VSA_EState4", "VSA_EState5", "VSA_EState6", "VSA_EState7", "VSA_EState8", "VSA_EState9", "VSA_EState10"]


# 파일경로 직접 설정
input_file = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\NPASS_NPs.xlsx"
output_file = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\output\\RDKit_output.xlsx"


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


#Descriptor 계산
rdkit_desc = []
for bi in range(0, len(mol_list)):

    print("count:", bi)
    print(smiles_list[bi])

    # 구조 생성문제 확인하기 위한 작업, 구조생성이 안된 경우 아무런 데이터도 포함되어 있지 않음
    type_check = str(mol_list[bi])

    if (type_check == 'None'):
        # Smile 정보로 mol 생성이 불가능한 경우 descriptor 계산이 안되므로 모두 0으로 처리 (에러방지)
        desc = []
        for bii in range(0, len(rdkit_desc_list)):
            desc.append(0)

        final_desc = []  # 각 분자별 descriptor를 최종정리 (리스트 내부의 리스트를 풀어서 하나의 리스트로 작성)
        for item in desc:
            if isinstance(item, list):
                final_desc.extend(item)
            else:
                final_desc.append(item)

        rdkit_desc.append(final_desc)

    else:
        desc = []  # 각 분자별 descriptor를 추가하기 위한 리스트 선언

        desc.append(Descriptors.BalabanJ(mol_list[bi]))  # BalabanJ
        desc.append(Descriptors.BertzCT(mol_list[bi]))  # BertzCT
        desc.append(Descriptors.Ipc(mol_list[bi]))  # Ipc
        desc.append(Descriptors.HallKierAlpha(mol_list[bi]))  # HallKierAlpha
        desc.append(Descriptors.Kappa1(mol_list[bi]))  # Kappa1
        desc.append(Descriptors.Kappa2(mol_list[bi]))  # Kappa2
        desc.append(Descriptors.Kappa3(mol_list[bi]))  # Kappa3
        desc.append(Descriptors.rdMolDescriptors.CalcPhi(mol_list[bi]))  # Phi
        desc.append(Descriptors.Chi0(mol_list[bi]))  # Chi0
        desc.append(Descriptors.Chi1(mol_list[bi]))  # Chi1
        desc.append(Descriptors.Chi0n(mol_list[bi]))  # Chi0n
        desc.append(Descriptors.Chi1n(mol_list[bi]))  # Chi1n
        desc.append(Descriptors.Chi2n(mol_list[bi]))  # Chi2n
        desc.append(Descriptors.Chi3n(mol_list[bi]))  # Chi3n
        desc.append(Descriptors.Chi4n(mol_list[bi]))  # Chi4n
        desc.append(Descriptors.Chi0v(mol_list[bi]))  # Chi0v
        desc.append(Descriptors.Chi1v(mol_list[bi]))  # Chi1v
        desc.append(Descriptors.Chi2v(mol_list[bi]))  # Chi2v
        desc.append(Descriptors.Chi3v(mol_list[bi]))  # Chi3v
        desc.append(Descriptors.Chi4v(mol_list[bi]))  # Chi4v
        desc.append(Descriptors.MolLogP(mol_list[bi]))  # MolLogP
        desc.append(Descriptors.MolMR(mol_list[bi]))  # MolMR
        desc.append(Descriptors.MolWt(mol_list[bi]))  # MolWt
        desc.append(Descriptors.ExactMolWt(mol_list[bi]))  # ExactMolWt
        desc.append(Descriptors.HeavyAtomCount(mol_list[bi]))  # HeavyAtomCount
        desc.append(Descriptors.HeavyAtomMolWt(mol_list[bi]))  # HeavyAtomMolWt
        desc.append(Descriptors.NHOHCount(mol_list[bi]))  # NHOHCount
        desc.append(Descriptors.NOCount(mol_list[bi]))  # NOCount
        desc.append(Descriptors.NumHAcceptors(mol_list[bi]))  # NumHAcceptors
        desc.append(Descriptors.NumHDonors(mol_list[bi]))  # NumHDonors
        desc.append(Descriptors.NumHeteroatoms(mol_list[bi]))  # NumHeteroatoms
        desc.append(Descriptors.NumRotatableBonds(mol_list[bi]))  # NumRotatableBonds
        desc.append(Descriptors.NumValenceElectrons(mol_list[bi]))  # NumValenceElectrons
        desc.append(Descriptors.rdMolDescriptors.CalcNumAmideBonds(mol_list[bi]))  # NumAmideBonds
        desc.append(Descriptors.rdMolDescriptors.CalcNumAromaticRings(mol_list[bi]))  # NumAromaticRings
        desc.append(Descriptors.rdMolDescriptors.CalcNumSaturatedRings(mol_list[bi]))  # NumSaturatedRings
        desc.append(Descriptors.rdMolDescriptors.CalcNumAliphaticRings(mol_list[bi]))  # NumAliphaticRings
        desc.append(Descriptors.rdMolDescriptors.CalcNumAromaticCarbocycles(mol_list[bi]))  # NumAromaticCarbocycles
        desc.append(Descriptors.rdMolDescriptors.CalcNumAromaticHeterocycles(mol_list[bi]))  # NumAromaticHeterocycles
        desc.append(Descriptors.rdMolDescriptors.CalcNumSaturatedCarbocycles(mol_list[bi]))  # NumSaturatedCarbocycles
        desc.append(Descriptors.rdMolDescriptors.CalcNumSaturatedHeterocycles(mol_list[bi]))  # NumSaturatedHeterocycles
        desc.append(Descriptors.rdMolDescriptors.CalcNumAliphaticCarbocycles(mol_list[bi]))  # NumAliphaticCarbocycles
        desc.append(Descriptors.rdMolDescriptors.CalcNumAliphaticHeterocycles(mol_list[bi]))  # NumAliphaticHeterocycles
        desc.append(Descriptors.RingCount(mol_list[bi]))  # RingCount
        desc.append(Descriptors.FractionCSP3(mol_list[bi]))  # FractionCSP3
        desc.append(Descriptors.rdMolDescriptors.CalcNumSpiroAtoms(mol_list[bi]))  # NumSpiroAtoms
        desc.append(Descriptors.rdMolDescriptors.CalcNumBridgeheadAtoms(mol_list[bi]))  # NumBridgeheadAtoms
        desc.append(Descriptors.TPSA(mol_list[bi]))  # TPSA
        desc.append(Descriptors.LabuteASA(mol_list[bi]))  # LabuteASA
        desc.append(Descriptors.rdMolDescriptors.PEOE_VSA_(mol_list[bi]))  # PEOE_VSA
        desc.append(Descriptors.rdMolDescriptors.SMR_VSA_(mol_list[bi]))  # SMR_VSA
        desc.append(Descriptors.rdMolDescriptors.SlogP_VSA_(mol_list[bi]))  # SlogP_VSA
        desc.append(Descriptors.EState_VSA1(mol_list[bi]))  # EState_VSA 1 ~ 11
        desc.append(Descriptors.EState_VSA2(mol_list[bi]))  # EState_VSA 1 ~ 11
        desc.append(Descriptors.EState_VSA3(mol_list[bi]))  # EState_VSA 1 ~ 11
        desc.append(Descriptors.EState_VSA4(mol_list[bi]))  # EState_VSA 1 ~ 11
        desc.append(Descriptors.EState_VSA5(mol_list[bi]))  # EState_VSA 1 ~ 11
        desc.append(Descriptors.EState_VSA6(mol_list[bi]))  # EState_VSA 1 ~ 11
        desc.append(Descriptors.EState_VSA7(mol_list[bi]))  # EState_VSA 1 ~ 11
        desc.append(Descriptors.EState_VSA8(mol_list[bi]))  # EState_VSA 1 ~ 11
        desc.append(Descriptors.EState_VSA9(mol_list[bi]))  # EState_VSA 1 ~ 11
        desc.append(Descriptors.EState_VSA10(mol_list[bi]))  # EState_VSA 1 ~ 11
        desc.append(Descriptors.EState_VSA11(mol_list[bi]))  # EState_VSA 1 ~ 11
        desc.append(Descriptors.VSA_EState1(mol_list[bi]))  # VSA_EState 1 ~ 10
        desc.append(Descriptors.VSA_EState2(mol_list[bi]))  # VSA_EState 1 ~ 10
        desc.append(Descriptors.VSA_EState3(mol_list[bi]))  # VSA_EState 1 ~ 10
        desc.append(Descriptors.VSA_EState4(mol_list[bi]))  # VSA_EState 1 ~ 10
        desc.append(Descriptors.VSA_EState5(mol_list[bi]))  # VSA_EState 1 ~ 10
        desc.append(Descriptors.VSA_EState6(mol_list[bi]))  # VSA_EState 1 ~ 10
        desc.append(Descriptors.VSA_EState7(mol_list[bi]))  # VSA_EState 1 ~ 10
        desc.append(Descriptors.VSA_EState8(mol_list[bi]))  # VSA_EState 1 ~ 10
        desc.append(Descriptors.VSA_EState9(mol_list[bi]))  # VSA_EState 1 ~ 10
        desc.append(Descriptors.VSA_EState10(mol_list[bi]))  # VSA_EState 1 ~ 10

        final_desc = []  # 각 분자별 descriptor를 최종정리 (리스트 내부의 리스트를 풀어서 하나의 리스트로 작성)

        for item in desc:
            if isinstance(item, list):
                final_desc.extend(item)
            else:
                final_desc.append(item)

        rdkit_desc.append(final_desc)

result = np.array(rdkit_desc)
result_df = pd.DataFrame(data=result, columns=rdkit_desc_list)
# result df에 CID, label column 추가
result_df.insert(loc=0, column="CAS", value=cas_list)
result_df.insert(loc=1, column="PubChem_CID", value=cid_list)
result_df.to_excel(output_file, index=False)