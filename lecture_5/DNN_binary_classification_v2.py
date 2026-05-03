# -*- coding: utf-8 -*-
# Author: Myungwon Seo

import os
import random
import numpy as np
import pandas as pd
import tensorflow as tf
import tensorflow.keras as keras
import matplotlib.pyplot as plt

from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    confusion_matrix,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    average_precision_score,
    balanced_accuracy_score,
    matthews_corrcoef
)
from imblearn.over_sampling import SMOTE
from tensorflow.keras.callbacks import EarlyStopping
from collections import Counter

# =========================================================
# 1. Seed 설정
# =========================================================
seed = 2023
random.seed(seed)
np.random.seed(seed)
tf.random.set_seed(seed)
os.environ["PYTHONHASHSEED"] = str(seed)

# =========================================================
# 2. 파일 경로
# =========================================================
input_file = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\NPASS_NPT204_RDKit_activity_output.xlsx"
output_file = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\lecture_5_output\\NPASS_NPT204_RDKit_activity_output_DNN.xlsx"
plot_dir = "C:\\Users\\mwseo\\PycharmProjects\\LAIDD_SMW\\dataset\\lecture_5_output\\"

os.makedirs(plot_dir, exist_ok=True)

# =========================================================
# 3. 데이터 불러오기
# =========================================================
df = pd.read_excel(input_file)

X = df.iloc[:, 5:111].values
y = df["label_classification"].values.astype(int)

# =========================================================
# 4. Train / Test split
#    중요: stratify 적용
# =========================================================
train_input, test_input, train_label, test_label = train_test_split(
    X,
    y,
    test_size=0.2,
    random_state=seed,
    stratify=y
)

# =========================================================
# 5. Normalization
#    중요: train에만 fit, test에는 transform
# =========================================================
scaler = MinMaxScaler()
train_input = scaler.fit_transform(train_input)
test_input = scaler.transform(test_input)


print("Before SMOTE:", Counter(train_label))
n0 = Counter(train_label)[0]
n1 = Counter(train_label)[1]

minority = min(n0, n1)
majority = max(n0, n1)

current_ratio = minority / majority
print("Current minority/majority ratio:", current_ratio)

# =========================================================
# 6. SMOTE
#    중요: train set에만 적용
# =========================================================
#smote = SMOTE(random_state=seed, sampling_strategy="auto")
smote = SMOTE(random_state=seed, sampling_strategy=0.5)
train_input, train_label = smote.fit_resample(train_input, train_label)

# float 변환
train_input = train_input.astype("float32")
test_input = test_input.astype("float32")
train_label = train_label.astype("float32")
test_label = test_label.astype("float32")

# =========================================================
# 7. 하이퍼파라미터
# =========================================================
learning_rate = 0.0001
n_epoch = 200
n_batch = 16
n_features = train_input.shape[1]

# =========================================================
# 8. Dataset 구성
#    repeat(), steps_per_epoch, validation_steps 사용하지 않음
# =========================================================
train_dataset = (
    tf.data.Dataset.from_tensor_slices((train_input, train_label))
    .shuffle(buffer_size=len(train_input), seed=seed)
    .batch(n_batch)
    .prefetch(tf.data.AUTOTUNE)
)

test_dataset = (
    tf.data.Dataset.from_tensor_slices((test_input, test_label))
    .batch(n_batch)
    .prefetch(tf.data.AUTOTUNE)
)

# =========================================================
# 9. DNN 모델 생성
#    이진분류: 출력 노드 1개 + sigmoid
# =========================================================
def create_model():
    model = keras.Sequential()

    model.add(keras.layers.Input(shape=(n_features,)))

    model.add(keras.layers.Dense(128, activation="relu"))
    model.add(keras.layers.Dropout(0.1))

    model.add(keras.layers.Dense(128, activation="relu"))
    model.add(keras.layers.Dropout(0.1))

    model.add(keras.layers.Dense(128, activation="relu"))
    model.add(keras.layers.Dropout(0.1))

    model.add(keras.layers.Dense(1, activation="sigmoid"))

    return model

model = create_model()

model.compile(
    optimizer=tf.keras.optimizers.Adam(learning_rate=learning_rate),
    loss="binary_crossentropy",
    metrics=[
        "accuracy",
        tf.keras.metrics.AUC(name="auc")
    ]
)

model.summary()

# =========================================================
# 10. EarlyStopping
# =========================================================
early_stop = EarlyStopping(
    monitor="val_loss",
    patience=30,
    restore_best_weights=True,
    verbose=1
)

# =========================================================
# 11. 모델 학습
# =========================================================
history = model.fit(
    train_dataset,
    epochs=n_epoch,
    validation_data=test_dataset,
    callbacks=[early_stop]
)

# =========================================================
# 12. 학습 곡선 저장
# =========================================================

# Accuracy
plt.figure()
plt.plot(history.history["accuracy"], label="Train Accuracy")
plt.plot(history.history["val_accuracy"], label="Validation Accuracy")
plt.xlabel("Epoch")
plt.ylabel("Accuracy")
plt.title("Train / Validation Accuracy")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(plot_dir, "accuracy_curve.png"), dpi=300, bbox_inches="tight")
plt.show()

# AUC
plt.figure()
plt.plot(history.history["auc"], label="Train AUC")
plt.plot(history.history["val_auc"], label="Validation AUC")
plt.xlabel("Epoch")
plt.ylabel("AUC")
plt.title("Train / Validation AUC")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(plot_dir, "auc_curve.png"), dpi=300, bbox_inches="tight")
plt.show()

# Loss
plt.figure()
plt.plot(history.history["loss"], label="Train Loss")
plt.plot(history.history["val_loss"], label="Validation Loss")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("Train / Validation Loss")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(plot_dir, "loss_curve.png"), dpi=300, bbox_inches="tight")
plt.show()

# Best epoch 표시 Loss 그래프
best_epoch = np.argmin(history.history["val_loss"]) + 1
best_val_loss = np.min(history.history["val_loss"])

plt.figure()
plt.plot(history.history["loss"], label="Train Loss")
plt.plot(history.history["val_loss"], label="Validation Loss")
plt.axvline(best_epoch, linestyle="--", label=f"Best Epoch: {best_epoch}")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("Train / Validation Loss with Best Epoch")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(plot_dir, "loss_curve_with_best_epoch.png"), dpi=300, bbox_inches="tight")
plt.show()

print(f"Best epoch based on val_loss: {best_epoch}")
print(f"Best val_loss: {best_val_loss:.4f}")

# =========================================================
# 13. 예측 및 성능 평가
# =========================================================
threshold = 0.5

train_prob = model.predict(train_input).ravel()
test_prob = model.predict(test_input).ravel()

train_pred = (train_prob >= threshold).astype(int)
test_pred = (test_prob >= threshold).astype(int)

# Confusion matrix
test_cm = confusion_matrix(test_label.astype(int), test_pred)
tn, fp, fn, tp = test_cm.ravel()

# Train metrics
train_acc = accuracy_score(train_label, train_pred)
train_auc = roc_auc_score(train_label, train_prob)
train_prauc = average_precision_score(train_label, train_prob)

# Test metrics
test_acc = accuracy_score(test_label, test_pred)
test_precision = precision_score(test_label, test_pred, zero_division=0)
test_auc = roc_auc_score(test_label, test_prob)
test_prauc = average_precision_score(test_label, test_prob)
test_sen = recall_score(test_label, test_pred, zero_division=0)
test_spe = tn / (tn + fp) if (tn + fp) > 0 else 0
test_ba = balanced_accuracy_score(test_label, test_pred)
test_f1 = f1_score(test_label, test_pred, zero_division=0)
test_mcc = matthews_corrcoef(test_label, test_pred)

# =========================================================
# 14. 결과 출력
# =========================================================
print("\n=== Train Performance ===")
print("Train ACC:", round(train_acc, 4))
print("Train AUC:", round(train_auc, 4))
print("Train PRAUC:", round(train_prauc, 4))

print("\n=== Test Performance ===")
print("Test ACC:", round(test_acc, 4))
print("Test Precision:", round(test_precision, 4))
print("Test AUC:", round(test_auc, 4))
print("Test PRAUC:", round(test_prauc, 4))
print("Test Sensitivity:", round(test_sen, 4))
print("Test Specificity:", round(test_spe, 4))
print("Test Balanced Accuracy:", round(test_ba, 4))
print("Test F1:", round(test_f1, 4))
print("Test MCC:", round(test_mcc, 4))

print("\nConfusion Matrix")
print(test_cm)

# =========================================================
# 15. 결과 저장
# =========================================================
final_df = pd.DataFrame({
    "Random Seed": [seed],
    "Threshold": [threshold],
    "Best Epoch": [best_epoch],
    "Best Val Loss": [best_val_loss],
    "Train ACC": [train_acc],
    "Train AUC": [train_auc],
    "Train PRAUC": [train_prauc],
    "Test ACC": [test_acc],
    "Test Precision": [test_precision],
    "Test AUC": [test_auc],
    "Test PRAUC": [test_prauc],
    "Test Sensitivity": [test_sen],
    "Test Specificity": [test_spe],
    "Test Balanced Accuracy": [test_ba],
    "Test F1": [test_f1],
    "Test MCC": [test_mcc],
    "TN": [tn],
    "FP": [fp],
    "FN": [fn],
    "TP": [tp]
})

final_df.to_excel(output_file, index=False)

print("\n결과 저장 완료:")
print(output_file)