# -*- coding: utf-8 -*-
# Author: Myungwon Seo
# Date: 2023-10-20
# E-mail: mwseo@krict.re.kr; seomyungwon@gmail.com

import os
import pandas as pd
import pubchempy as pcp

# 파일경로 설정
current_directory = os.path.dirname(__file__)
input_file = os.path.join(current_directory, "..", "dataset", "NPASS_NPs.xlsx")

# 파일경로 직접 설정
#input_file = "userpath\\NPASS_NPs.xlsx"
input_file = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\NPASS_NPs.xlsx"


#파일을 dataframe 형태로 불러오기
df = pd.read_excel(input_file)

#파일속성 확인
df.columns

#CID 값 추출
cid = df['PubChem_CID']
cid_single = cid[0]

#분자 불러오기
c = pcp.Compound.from_cid(str(cid_single)) #CID 값을 문자열로 입력해야 함
c = pcp.get_compounds('Eupatilin', 'name') #이름으로 검색
c = pcp.get_compounds('COC1=C(C=C(C=C1)C2=CC(=O)C3=C(O2)C=C(C(=C3O)OC)O)OC', 'smiles') #SMIELS로 검색
c = pcp.get_compounds('Eupatilin', 'name', record_type='3d') # 3D 구조 불러옴

print(df)

#분자 특성 불러오기
c_mf = c.molecular_formula # molecular formula
c_mw = c.molecular_weight #molecular weight
c_smiles = c.isomeric_smiles #SMILES
c_xlogp = c.xlogp #logP
c_name = c.iupac_name # iupac name
c_synnonyms = c.synonyms # synonyms

pcp.get_properties('AtomStereoCount', '6623')

#사전형태로 데이터 불러오기 #get_properties('Property', 'identifier(CID)')
p = pcp.get_properties('IsomericSMILES', 'CC', 'smiles', searchtype='superstructure')
p = pcp.get_properties('XLogP', '6623')



''' Pubchem properties
MolecularFormula, MolecularWeight, CanonicalSMILES, IsomericSMILES, 
InChI, InChIKey, IUPACName, XLogP, ExactMass, MonoisotopicMass, TPSA, 
Complexity, Charge, HBondDonorCount, HBondAcceptorCount, RotatableBondCount, 
HeavyAtomCount, IsotopeAtomCount, AtomStereoCount, DefinedAtomStereoCount, 
UndefinedAtomStereoCount, BondStereoCount, DefinedBondStereoCount, UndefinedBondStereoCount, 
CovalentUnitCount, Volume3D, XStericQuadrupole3D, YStericQuadrupole3D, ZStericQuadrupole3D, 
FeatureCount3D, FeatureAcceptorCount3D, FeatureDonorCount3D, FeatureAnionCount3D, FeatureCationCount3D, 
FeatureRingCount3D, FeatureHydrophobeCount3D, ConformerModelRMSD3D, EffectiveRotorCount3D, ConformerCount3D.
'''


