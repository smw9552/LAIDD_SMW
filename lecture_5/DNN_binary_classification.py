# -*- coding: utf-8 -*-
# Author: Myungwon Seo
# Date: 2023-10-20
# E-mail: mwseo@krict.re.kr; seomyungwon@gmail.com

import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from imblearn.over_sampling import SMOTE
import tensorflow.keras as keras
import tensorflow as tf
from sklearn.metrics import confusion_matrix
import random
from sklearn.model_selection import StratifiedShuffleSplit

# 파일경로 설정
#current_directory = os.path.dirname(__file__)
#input_file = os.path.join(current_directory, "..", "dataset", "output", "NPASS_NPT204_RDKit_activity_output.xlsx")
#output_file = os.path.join(current_directory, "..", "dataset", "output", "NPASS_NPT204_RDKit_activity_output_DNN.xlsx")

#seed 설정(재현성을 위한 작업)
seed = 2023
random.seed(seed)
os.environ["PYTHONHASHSEED"] = str(seed)
tf.random.set_seed(seed)

# 파일경로 직접 설정
input_file = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\NPASS_NPT204_RDKit_activity_output.xlsx"
output_file = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\NPASS_NPT204_RDKit_activity_output_DNN.xlsx"

#perforamnce list
train_accuracy = []
train_AUC = []
test_accuracy = []
test_AUC = []
test_sensitivity = []
test_specificity = []
test_balanced_accuracy = []


#파일을 dataframe 형태로 불러오기
df = pd.read_excel(input_file)
df_data = df.iloc[:,5:111] #df에서 전처리 필요한 부분만 추출
label = df['label_classification'] # 데이터 label을 리스트로 불러옴

# 정규화 (Normalization)
scaler = MinMaxScaler()
df_data_new = scaler.fit_transform(df_data)

# 딥러닝 훈련 데이터 입력을 위한 형태 변환
data = np.array(df_data_new)
target = np.array(label)


# 데이터 분할 (일반적인 방법), StratifiedShuffleSplit을 활용하면 label의 분포를 고려하여 선별 가능
train_input, test_input, train_label, test_label = train_test_split(data, target, test_size = 0.2, random_state=seed)


# 데이터 분할 (label의 분포를 고려)
'''
n_splits = 10
test_size = 0.2  # 20%
random_state = seed
splitter = StratifiedShuffleSplit(n_splits=n_splits, test_size=test_size, random_state=random_state)

for train_index, test_index in splitter.split(data, target):
    train_input, test_input = data[train_index], data[test_index]
    train_label, test_label = target[train_index], target[test_index]
'''

# 데이터 샘플링 (SMOTE)
smote = SMOTE(random_state=seed, sampling_strategy='auto')
# smote = SMOTE(random_state=seed, sampling_strategy=0.6, sampling_strategy = 'auto')
train_input_smote, train_label_smote = smote.fit_resample(train_input, train_label)

train_input = train_input_smote
train_label = train_label_smote

# 하이퍼파라미터
learning_rate = 0.0001 # 0.001 default: 0.001
n_epoch = 400 # 300
n_batch = 16  # 16
n_class = 2
n_train = train_input.shape[0]  # training 데이터 수
n_features = train_input.shape[1]  # feature(descriptor) 수
n_test = test_input.shape[0]  # test 데이터 수

# one-hot encoding
train_label = keras.utils.to_categorical(train_label, n_class)  # label encoding of training data
test_label = keras.utils.to_categorical(test_label, n_class)  # label encoding of test data

# Dataset 구성(desc, label을 묶어서 제공)
train_dataset = tf.data.Dataset.from_tensor_slices((train_input, train_label)).shuffle(100000).batch(n_batch, drop_remainder=True).repeat()
test_dataset = tf.data.Dataset.from_tensor_slices((test_input, test_label)).batch(n_batch)


# 모델 생성
def create_model():
    model = keras.Sequential()

    # 입력층 및 hidden layer
    model.add(keras.layers.Dense(350, activation='relu', input_shape=(n_features,)))
    model.add(keras.layers.Dropout(0.5))
    model.add(keras.layers.Dense(350, activation='relu'))
    model.add(keras.layers.Dropout(0.5))
    model.add(keras.layers.Dense(350, activation='relu'))
    model.add(keras.layers.Dropout(0.5))
    model.add(keras.layers.Dense(350, activation='relu'))
    model.add(keras.layers.Dropout(0.5))

    # optimized
    # model.add(keras.layers.Dense(350, activation='relu', input_shape=(n_features,)))
    # model.add(keras.layers.Dropout(0.5))
    # model.add(keras.layers.Dense(350, activation='relu'))
    # model.add(keras.layers.Dropout(0.5))
    # model.add(keras.layers.Dense(350, activation='relu'))
    # model.add(keras.layers.Dropout(0.5))
    # model.add(keras.layers.Dense(350, activation='relu'))
    # model.add(keras.layers.Dropout(0.5))
    # model.add(keras.layers.Dense(350, activation='relu'))
    # model.add(keras.layers.Dropout(0.5))
    # model.add(keras.layers.Dense(350, activation='relu'))
    # model.add(keras.layers.Dropout(0.5))


    # 출력층 #
    # activation = softmax (multi class), activtion = sigmoid (binary class)
    model.add(keras.layers.Dense(n_class, activation='sigmoid'))

    return model

# 모델생성 & 컴파일
# loss 함수
# multi-class: categorical_crossentropy
# binary-class: binary_crossentropy
# (참고) metric: https://keras.io/api/metrics/

model = create_model()
#model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate), loss='binary_crossentropy', metrics=['accuracy', 'AUC'])

model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate),
              loss='binary_crossentropy',
              metrics=['accuracy', 'AUC', 'TruePositives', 'TrueNegatives', 'FalsePositives', 'FalseNegatives'])

steps_per_epoch = n_train // n_batch
validation_steps = int(np.ceil(n_test / n_batch))

# 학습
history = model.fit(train_dataset,
                    epochs=n_epoch,
                    steps_per_epoch=steps_per_epoch,
                    validation_data=test_dataset,
                    validation_steps=validation_steps)

# 예측모델 최종 80% 내부검증(성능)
train_acc = history.history['accuracy'][n_epoch - 1]
train_auc = history.history['auc'][n_epoch - 1]

# 예측모델 최종 20% 외부검증(성능)
model.evaluate(test_dataset)

# 예측모델 최종 20% 외부검증(예측 + one-hot encoding)
test_label_predict = model.predict(test_dataset)
test_label_pred = keras.utils.to_categorical(np.argmax(test_label_predict, axis=1), n_class)

#confusion matrix 계산을 위한 decoding
test_label_decoding = []
test_label_pred_decoding = []
for ci in range(0, len(test_label)): #test_label, test_label_pred 숫자가 동일함
    test_label_decoding.append(np.argmax(test_label[ci]))
    test_label_pred_decoding.append(np.argmax(test_label_pred[ci]))

test_label_decode = np.array(test_label_decoding)
test_label_pred_decode = np.array(test_label_pred_decoding)
test_cm = confusion_matrix(test_label_decode, test_label_pred_decode) #decoding 형태로 데이터를 추가해야 confusion matrix 계산 가능
#(참고) test_cm[0,0] = true negative, test_cm[0,1] = false_positive, test_cm[1,0] = false_negative, test_cm[1,1] = true_positive
# TN | FP
# FB | TP

# test 평가지표 계산
test_tp = test_cm[1,1]
test_tn = test_cm[0,0]
test_fp = test_cm[0,1]
test_fn = test_cm[1,0]

test_acc = (test_tp + test_tn) / (test_tp + test_tn + test_fp + test_fn)
test_auc = history.history['val_auc'][n_epoch - 1]
test_sen = test_tp / (test_tp + test_fn)
test_spe = test_tn / (test_tn + test_fp)
test_ba = (test_sen + test_spe) / 2

# 성능평가지표 출력
print("Train acc: ", round(train_acc, 3))
print("Train auc: ", round(train_auc, 3))

print("Test acc: ", test_acc)
print("Test auc: ", test_auc)
print("Test sen: ", test_sen)
print("Test spe: ", test_spe)
print("Test ba: ", test_ba)


#결과파일 작성
train_accuracy.append(train_acc)
train_AUC.append(train_auc)

test_accuracy.append(test_acc)
test_AUC.append(test_auc)
test_sensitivity.append(test_sen)
test_specificity.append(test_spe)
test_balanced_accuracy.append(test_ba)

final_df = pd.DataFrame({'Random Seed': seed,
                         'Internal ACC': train_accuracy,
                         'Internal AUC': train_AUC,
                         'External ACC': test_accuracy,
                         'External AUC': test_AUC,
                         'External Sensitivity': test_sensitivity,
                         'External Specificity': test_specificity,
                         'External BA': test_balanced_accuracy,
                         })

final_df.to_excel(output_file)