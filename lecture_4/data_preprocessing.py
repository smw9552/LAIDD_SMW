# -*- coding: utf-8 -*-
# Author: Myungwon Seo
# Date: 2023-10-20
# E-mail: mwseo@krict.re.kr; seomyungwon@gmail.com


import pandas as pd
import os
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from imblearn.under_sampling import RandomUnderSampler
from imblearn.over_sampling import SMOTE

# 파일경로 설정
current_directory = os.path.dirname(__file__)
input_file = os.path.join(current_directory, "..", "dataset", "output", "NPASS_NPT204_RDKit_activity_output.xlsx")
output_file_nor = os.path.join(current_directory, "..", "dataset", "output", "NPASS_NPT204_RDKit_activity_output_nor.xlsx")
output_file_std = os.path.join(current_directory, "..", "dataset", "output", "NPASS_NPT204_RDKit_activity_output_std.xlsx")

# 파일경로 직접 설정
input_file = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\NPASS_NPT204_RDKit_activity_output.xlsx"
output_file_nor = "/dataset/output/NPASS_NPT204_RDKit_activity_output_nor.xlsx"
output_file_std = "/dataset/output/NPASS_NPT204_RDKit_activity_output_std.xlsx"

#파일을 dataframe 형태로 불러오기
df = pd.read_excel(input_file)
df_data = df.iloc[:,5:111] #df에서 전처리 필요한 부분만 추출

# 정규화(Normalization)
scaler = MinMaxScaler()
normalization = scaler.fit_transform(df_data)
df_data_nor = pd.DataFrame(normalization, columns=df_data.columns)
df_data_nor.to_excel(output_file_nor, index=False) #데이터 확인

# 표준화(Standardization)
std_scaler = StandardScaler()
standardization = std_scaler.fit_transform(df_data)
df_data_std = pd.DataFrame(standardization, columns=df_data.columns)
df_data_std.to_excel(output_file_std, index=False) #데이터 확인


# 데이터 샘플링 (undersampling)
label = df['label_classification'] # 데이터 label을 리스트로 불러옴

data = np.array(df_data)
target = np.array(label)

#total
print(str("total data:"), len(target))
#positive
print(str("positive data:"), sum(target))
#negative
print(str("negative data:"), len(target) - sum(target))

# positive 데이터 수를 기준으로 negative 개수를 지정 (1:1 비율로 수정)
sampler = RandomUnderSampler(sampling_strategy={0: len(np.where(target == 1)[0]), 1: len(np.where(target == 1)[0])})
data_downsampled, target_downsampled = sampler.fit_resample(data, target)

#total
print(str("total data:"), len(target_downsampled))
#positive
print(str("positive data:"), sum(target_downsampled))
#negative
print(str("negative data:"), len(target_downsampled) - sum(target_downsampled))

# 데이터 샘플링 (oversampling)
seed = 2023
smote = SMOTE(random_state=seed, sampling_strategy='auto') #sampling_stragety는 smapling 비율 설정
data_smote, target_smote = smote.fit_resample(data, target)

#total
print(str("total data:"), len(target))
#positive
print(str("positive data:"), sum(target))
#negative
print(str("negative data:"), len(target) - sum(target))

#total
print(str("total data:"), len(target_smote))
#positive
print(str("positive data:"), sum(target_smote))
#negative
print(str("negative data:"), len(target_smote) - sum(target_smote))