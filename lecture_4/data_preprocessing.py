import pandas as pd
import os


# 파일경로 설정
current_directory = os.path.dirname(__file__)
input_file = os.path.join(current_directory, "..", "dataset", "output", "RDKit_output.xlsx")


# 파일경로 직접 설정
#input_file = "userpath\\NPASS_NPs.xlsx"
input_file = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\output\\RDKit_output.xlsx"
output_file = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\output\\Mordred_output.xlsx"

#파일을 dataframe 형태로 불러오기
df = pd.read_excel(input_file)

#파일속성 확인
df.columns