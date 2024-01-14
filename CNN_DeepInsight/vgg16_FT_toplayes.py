import numpy as np
import tensorflow as tf
from tensorflow import keras
import random


import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import os
from scipy.stats import zscore

from tensorflow.keras.preprocessing.image import array_to_img, img_to_array, load_img
from tensorflow.keras.applications import VGG16, imagenet_utils
# from keras.callbacks import EarlyStopping
from keras.callbacks import EarlyStopping, ModelCheckpoint

# model = VGG16(weights='imagenet')
# model.summary()

met_control = pd.read_csv("../../Met/met_rr_control_named.csv")
print(len(met_control))


labels = zscore(np.array(met_control.glycerol)) #normalization
# labels = np.array(met_control.glycerol)
# print(labels)
plt.hist(labels)
lines = np.array(met_control.NSFTV_ID)
# print(lines)

dff = pd.DataFrame({
        "labels": labels,
        "NSFTV_ID": lines
    })

dff.head()

dataset_path = "../../SNPimg/dataset/"
chr_path = os.listdir("../../SNPimg/dataset/NSFTV_1/")
base_df = dff.assign(dataset_path=dataset_path)
base_df.head()

image_pathss = [base_df["dataset_path"] + base_df["NSFTV_ID"] + "/" + chr_path[i] for i in range(12)]  ####### here we can include multiple chromosome setting.
# print(image_paths)
print(np.array(image_pathss).shape)
# print(image_paths)
print(np.array(image_pathss).shape)

image_total = []

# loop over each batch

for image_paths in image_pathss:
  # print(image_paths)
  batch_size = 162
  # image_features_chr = []
  image_labels = []
  # image_paths = image_pathss[0]
  for i in range(0, len(image_paths)//batch_size):
    # print(i)
    batch_paths = image_paths[i:i + batch_size]
    batch_labels = dff["labels"][i:i + batch_size]
    batch_images = []
    for image_path in batch_paths:
      image = load_img(image_path, target_size = (224, 224)) #image is still in PIL format. Need to convert into np.array
      image = img_to_array(image)
      # We expand the dimensions and then subtract the mean RGB pixel intensity of ImageNet
      image = np.expand_dims(image, axis=0) # add one dim in axis=0
      image = imagenet_utils.preprocess_input(image) #normalization based on ImageNet.
      batch_images.append(image)
    
    batch_images = np.vstack(batch_images) #looks like no need to expand in the before lines.
    # features = model.predict(batch_images, batch_size = batch_size)
    # print(features.shape)
    # features = np.reshape(features,(-1, features.shape[1]*features.shape[2]*features.shape[3]))
    # print(features.shape)
    # # store our features and corresponding labels
    # image_features_chr.append(features)
    image_labels.append(batch_labels)

  image_total.append(batch_images)

# image_features.append(np.array(image_features_chr))  

i = 1
i = str(i)
cvv = pd.read_csv(f'../../Met/CrossValidation/cv_{i}/met_cv_{i}.csv')
train_index = cvv.query("Treatment=='Control' and set == 'train'").index
test_index = cvv.query("Treatment=='Control' and set == 'test'").index

train_labels = batch_labels[train_index]
test_labels = batch_labels[test_index]

Train_image = []
Test_image = []
for chrr in range(0,12):
    train_images = image_total[chrr][train_index,:,:,:]
    test_images = image_total[chrr][test_index,:,:,:]
    
    Train_image.append(train_images)
    Test_image.append(test_images)

size = (150, 150)

def newdataset_ge(inputs, targets, size):
    inputs = tf.convert_to_tensor(inputs, dtype=tf.float32)  # Adjust dtype as needed
    targets = tf.convert_to_tensor(targets, dtype=tf.float32)  # Adjust dtype as needed

    # Create a dataset using from_tensor_slices
    new_dataset = tf.data.Dataset.from_tensor_slices(({"inputs1": tf.image.resize(inputs[0], size=size), 
                                                       "inputs2": tf.image.resize(inputs[1], size=size), 
                                                       "inputs3":tf.image.resize(inputs[2], size=size),
                                                       "inputs4": tf.image.resize(inputs[3], size=size), 
                                                       "inputs5": tf.image.resize(inputs[4], size=size), 
                                                       "inputs6":tf.image.resize(inputs[5], size=size),
                                                       "inputs7": tf.image.resize(inputs[6], size=size), 
                                                       "inputs8": tf.image.resize(inputs[7], size=size), 
                                                       "inputs9":tf.image.resize(inputs[8], size=size),
                                                       "inputs10": tf.image.resize(inputs[9], size=size), 
                                                       "inputs11": tf.image.resize(inputs[10], size=size), 
                                                       "inputs12":tf.image.resize(inputs[11], size=size)}, 
                                                       targets))
    return new_dataset

train_ds = newdataset_ge(Train_image, train_labels, size)
test_ds = newdataset_ge(Test_image, test_labels, size)


batch_size = 8

train_ds = train_ds.cache().batch(batch_size).prefetch(buffer_size=10)
# validation_ds = validation_ds.cache().batch(batch_size).prefetch(buffer_size=10)
test_ds = test_ds.cache().batch(batch_size).prefetch(buffer_size=10)

base_model_vgg = VGG16(
    weights="imagenet",  # Load weights pre-trained on ImageNet.
    input_shape=(150, 150, 3),
    include_top=False,
)  # Do not include the ImageNet classifier at the top.

# Freeze the base_model

base_model_vgg.trainable = False


# Create new model on top
inputs1 = keras.Input(shape=(150, 150, 3), name="inputs1")
inputs2 = keras.Input(shape=(150, 150, 3), name="inputs2")
inputs3 = keras.Input(shape=(150, 150, 3), name="inputs3")
inputs4 = keras.Input(shape=(150, 150, 3), name="inputs4")
inputs5 = keras.Input(shape=(150, 150, 3), name="inputs5")
inputs6 = keras.Input(shape=(150, 150, 3), name="inputs6")
inputs7 = keras.Input(shape=(150, 150, 3), name="inputs7")
inputs8 = keras.Input(shape=(150, 150, 3), name="inputs8")
inputs9 = keras.Input(shape=(150, 150, 3), name="inputs9")
inputs10 = keras.Input(shape=(150, 150, 3), name="inputs10")
inputs11 = keras.Input(shape=(150, 150, 3), name="inputs11")
inputs12 = keras.Input(shape=(150, 150, 3), name="inputs12")

def multi_branch_build(inputs, base_model):

    # inputs = keras.Input(shape=(150, 150, 3))
    # x = data_augmentation(inputs)  # Apply random data augmentation

    x = inputs
    # Pre-trained Xception weights requires that input be scaled
    # from (0, 255) to a range of (-1., +1.), the rescaling layer
    # outputs: `(inputs * scale) + offset`
    scale_layer = keras.layers.Rescaling(scale=1 / 127.5, offset=-1)
    x = scale_layer(x)

    # The base model contains batchnorm layers. We want to keep them in inference mode
    # when we unfreeze the base model for fine-tuning, so we make sure that the
    # base_model is running in inference mode here.
    x = base_model(x, training=False)
    # x=keras.layers.GlobalAveragePooling2D()(x)
  # Regularize with dropout
    return x

x1 = multi_branch_build(inputs = inputs1, base_model=base_model_vgg)
x2 = multi_branch_build(inputs = inputs2, base_model=base_model_vgg)
x3 = multi_branch_build(inputs = inputs3, base_model=base_model_vgg)
x4 = multi_branch_build(inputs = inputs4, base_model=base_model_vgg)
x5 = multi_branch_build(inputs = inputs5, base_model=base_model_vgg)
x6 = multi_branch_build(inputs = inputs6, base_model=base_model_vgg)
x7 = multi_branch_build(inputs = inputs7, base_model=base_model_vgg)
x8 = multi_branch_build(inputs = inputs8, base_model=base_model_vgg)
x9 = multi_branch_build(inputs = inputs9, base_model=base_model_vgg)
x10 = multi_branch_build(inputs = inputs10, base_model=base_model_vgg)
x11 = multi_branch_build(inputs = inputs11, base_model=base_model_vgg)
x12 = multi_branch_build(inputs = inputs12, base_model=base_model_vgg)

y = keras.layers.concatenate([x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12])
y = keras.layers.GlobalAveragePooling2D()(y)
y = keras.layers.Dropout(0.2)(y)
y = keras.layers.Dense(256, activation='relu')(y)
y = keras.layers.Dropout(0.2)(y)
y = keras.layers.Dense(128, activation='relu')(y)
outputs = keras.layers.Dense(1)(y)
model = keras.Model([inputs1, inputs2, inputs3, inputs4, inputs5, inputs6, inputs7, inputs8, inputs9, inputs10, inputs11, inputs12], outputs)

model.summary()

seed = 42
random.seed(seed)
np.random.seed(seed)
tf.random.set_seed(seed)
keras.utils.set_random_seed(seed)


model.compile(
    optimizer=tf.keras.optimizers.legacy.Adam(learning_rate=1e-3),
    loss=keras.losses.MeanSquaredError(),
    metrics=[keras.metrics.MeanSquaredError()],
    # batch_size=32
)

es = EarlyStopping(monitor='val_loss', patience=100)
epochs = 500

model.fit(train_ds.take(130).cache(), 
        epochs=epochs, 
        validation_data=test_ds,
        callbacks = [es])

plt.plot(model.history.history['loss'], label='loss')
plt.plot(model.history.history['val_loss'], label = 'val_loss')

plt.xlabel('Epoch')
plt.ylabel('LOSS')
plt.legend(loc='upper right')
plt.show()

y_pred = model.predict(test_ds)

from scipy.stats import pearsonr
correlation, _ = pearsonr(np.array(y_pred.reshape(32)), np.array(test_labels))
print(correlation)