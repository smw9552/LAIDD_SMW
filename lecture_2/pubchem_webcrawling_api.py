# -*- coding: utf-8 -*-
# Author: Myungwon Seo
# Date: 2023-10-20
# E-mail: mwseo@krict.re.kr; seomyungwon@gmail.com

import os
import pandas as pd
from urllib.error import HTTPError
from urllib.request import urlopen
from bs4 import BeautifulSoup
import requests


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

#CAS 정보 추출하기
cas = df['CAS']


## PubChem WebCrawling ##

#CAS 정보만 있을 때 PubChem에서 CID 정보 가지고 오는 코드
pubchem_cid = []
for cas_num in cas:

    try:
        label = True

        print(str("Input CAS number: ") + str(cas_num).strip())
        pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/compound/" + str(cas_num).strip()
        cid_info = urlopen(pubchem_url, None, timeout=10000)

        while True:
            pubchem_line = cid_info.readline()
            if not pubchem_line: break

            if (str(pubchem_line).__contains__('''b'    <meta name="pubchem_uid_value"''')):
                pubchem_cid.append(str(pubchem_line).replace("b'", "").
                                    replace('"', '').
                                    replace("    <meta name=pubchem_uid_value content=", "").
                                    replace("'", "").
                                    replace(">", "").
                                    replace("\\n", "").strip())


                label = False
                break

        if(label):
            pubchem_cid.append("None")

    except HTTPError as e:
        print(e)
        pubchem_cid.append("None")



#CID 정보로 PubChem에서 property 정보 가지고 오는 코드 (예시 canonical smiles)
smiles = []
for cid_num in pubchem_cid:

    try:
        label = True

        print(str("Input CID number: ") + str(cid_num).strip())
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + str(cid_num).strip() + "/property/CanonicalSMILES/txt"

        ''' 
        MolecularFormula, MolecularWeight, CanonicalSMILES, IsomericSMILES, 
        InChI, InChIKey, IUPACName, XLogP, ExactMass, MonoisotopicMass, TPSA, 
        Complexity, Charge, HBondDonorCount, HBondAcceptorCount, RotatableBondCount, 
        HeavyAtomCount, IsotopeAtomCount, AtomStereoCount, DefinedAtomStereoCount, 
        UndefinedAtomStereoCount, BondStereoCount, DefinedBondStereoCount, UndefinedBondStereoCount, 
        CovalentUnitCount, Volume3D, XStericQuadrupole3D, YStericQuadrupole3D, ZStericQuadrupole3D, 
        FeatureCount3D, FeatureAcceptorCount3D, FeatureDonorCount3D, FeatureAnionCount3D, FeatureCationCount3D, 
        FeatureRingCount3D, FeatureHydrophobeCount3D, ConformerModelRMSD3D, EffectiveRotorCount3D, ConformerCount3D.
        '''

        url_info = urlopen(url, None, timeout=100000)

        bsObject = BeautifulSoup(url_info, "html.parser")
        All_txt = bsObject.text

        smiles.append(All_txt.strip("\n"))

        label = False

        if(label):
            smiles.append("None")

    except HTTPError as e:
        print(e)
        smiles.append("None")



#CID 정보로 PubChem에서 assay 정보 가져오는 방법
#CID에 해당하는 assay list 불러오기 (CID가 포함된 Assay ID 전부 호출

aids = []
for cid_num in pubchem_cid:

    try:
        label = True

        print(str("Input CID number: ") + str(cid_num).strip())
        assay_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + str(cid_num).strip() + "/aids/txt"

        response = requests.get(assay_url, verify=False)


        aids.append(response.text.strip().split('\n'))

        label = False

        if(label):
            aids.append("None")

    except HTTPError as e:
        print(e)
        aids.append("None")



#하나의 assay ID로 실험 데이터 다운로드
target_aid = aids[0][0] #위에서 수집했던 aid 중 한개만 추출
assay_download_url = "https://pubchem.ncbi.nlm.nih.gov/assay/pcget.cgi?query=download&record_type=datatable&actvty=all&response_type=save&aid=" + str(target_aid)
download_response = requests.get(assay_download_url, verify=False)
filepath = "C:\\Users\\user\\Desktop\\test.csv"
with open(filepath, "wb") as file:
    file.write(download_response.content)



#다수의 assay ID로 실험데이터 추출
file_dir = "C:\\Users\\user\\Desktop\\"

for aid_list in aids:
    for aid_num in aid_list:

        try:
            assay_download_url = "https://pubchem.ncbi.nlm.nih.gov/assay/pcget.cgi?query=download&record_type=datatable&actvty=all&response_type=save&aid=" + str(aid_num)
            download_response = requests.get(assay_download_url, verify=False)
            download_file_info = str(file_dir) + str(aid_num) + str(".csv")
            with open(download_file_info, "wb") as file:
                file.write(download_response.content)

        except HTTPError as e:
            print(e)
            aids.append("None")