from tensorflow.keras.preprocessing.image import array_to_img, img_to_array, load_img
from tensorflow.keras.applications import VGG16, imagenet_utils
from tensorflow.keras.applications.inception_v3 import InceptionV3
from tensorflow.keras.applications.mobilenet_v2 import MobileNetV2
from tensorflow.keras.applications.densenet import DenseNet201
from tensorflow.keras.applications.nasnet import NASNetMobile



from tensorflow.keras.applications.resnet50 import ResNet50
from tensorflow.keras.applications.efficientnet import EfficientNetB7
from tensorflow import keras
import pandas as pd 
import numpy as np
import random
import sys
import os
from scipy.stats import zscore
from scipy.stats import pearsonr
# from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
# from sklearn.ensemble import RandomForestRegressor
import argparse

parser = argparse.ArgumentParser(description = 'Cross validation for DL4GP')
parser.add_argument("--cv", dest='CV', help = 'Cv start point', type=int)
parser.add_argument("--modell", dest='modell', help = 'DL model to extract features')

args = parser.parse_args()
cv_1 = args.CV
modell = args.modell

outfile = "log.out_" + modell + "_cv" + str(cv_1)


def image_extractor(model_name):
    if(model_name == "VGG16"):
        model = VGG16(weights="imagenet", include_top=False)
        # model.summary()

    if(model_name == "ResNet50"):
        model = ResNet50(weights="imagenet", include_top=False)
        # model.summary()

    if(model_name == "EfficientNetB7"):
        model = EfficientNetB7(weights="imagenet", include_top=False)
        # model.summary()
    if(model_name == "InceptionV3"):
        model = InceptionV3(weights='imagenet', include_top=False)

    if(model_name == "MobileNetV2"):
        model = MobileNetV2(weights='imagenet', include_top=False)

    if(model_name == "DenseNet201"):
        model = DenseNet201(weights='imagenet', include_top=False)

    if(model_name == "NASNetMobile"):
        model = NASNetMobile(weights='imagenet', include_top=False)

########################################################################
    met_control = pd.read_csv("../../Met/met_rr_control_named.csv")
    labels = zscore(np.array(met_control.glycerol)) #normalization
    lines = np.array(met_control.NSFTV_ID)

    dff = pd.DataFrame({
            "labels": labels,
            "NSFTV_ID": lines
        })

    dataset_path = "../../SNPimg/dataset/"
    chr_path = os.listdir("../../SNPimg/dataset/NSFTV_1/")
    base_df = dff.assign(dataset_path=dataset_path)

    image_pathss = [base_df["dataset_path"] + base_df["NSFTV_ID"] + "/" + chr_path[i] for i in range(12)]  ####### here we can include multiple chromosome setting.
    print(np.array(image_pathss).shape)

    image_features = []
    for image_paths in image_pathss:
        batch_size = 162
        image_features_chr = []
        image_labels = []

        for i in range(0, len(image_paths)//batch_size):
            batch_paths = image_paths[i:i + batch_size]
            batch_labels = dff["NSFTV_ID"][i:i + batch_size]
            batch_images = []
            for image_path in batch_paths:
                image = load_img(image_path, target_size = (244, 244)) #image is still in PIL format. Need to convert into np.array
                image = img_to_array(image)
                # We expand the dimensions and then subtract the mean RGB pixel intensity of ImageNet
                image = np.expand_dims(image, axis=0) # add one dim in axis=0
                image = imagenet_utils.preprocess_input(image) #normalization based on ImageNet.
                # print(np.min(image), np.max(image))
                image = image/127.5-1.0 #scale into range (-1,1)
                # print(np.min(image), np.max(image))
                batch_images.append(image)

            batch_images = np.vstack(batch_images) #looks like no need to expand in the before lines.
            features = model.predict(batch_images, batch_size = batch_size)
            print(features.shape)
            features = np.reshape(features,(-1, features.shape[1]*features.shape[2]*features.shape[3]))
            print(features.shape)
                # # store our features and corresponding labels
            image_features_chr.append(features)
            image_labels.append(batch_labels)
        
        image_features.append(np.array(image_features_chr))  
    
    image_featuress = np.concatenate(np.array(image_features), axis=2) # concatenate multiple chr together.


    # take our list of batches and reduce the dimernsion so that it's now a list 25088 features x 25000 rows (25000 x 1 for our labels)
    imageLabels_data =  [lb for label_batch in image_labels for lb in label_batch]
    imageFeatures_data = [feature for feature_batch in image_featuress for feature in feature_batch]

    # Convert to numpy arrays
    image_labels_data = np.array(imageLabels_data)
    image_features_data = np.array(imageFeatures_data)

    return image_features_data


def ll_reg(modell, image_features_data, i, alpha, outputFolder):

        i = str(i)
        cvv = pd.read_csv(f"../../Met/CrossValidation/cv_{i}/met_cv_{i}.csv")
        train_index = cvv.query("Treatment=='Control' and set == 'train'").index
        test_index = cvv.query("Treatment=='Control' and set == 'test'").index

        met_control = pd.read_csv("../../Met/met_rr_control_named.csv")
        met_stress = pd.read_csv("../../Met/met_rr_stress_named.csv")

        labels = zscore(np.hstack([met_control.iloc[:,2:], met_stress.iloc[:,2:]]))
        y = labels 
        X = image_features_data

        y_train = labels[train_index]
        y_test = labels[test_index]
        X_train = X[train_index,:]
        X_test = X[test_index,:]
        print(y_train.shape)
        print(X_train.shape)
        # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

        # lm = LinearRegression()
        # lm.fit(X_train, y_train)
        # y_pred = lm.predict(X_test)
        # regressor = RandomForestRegressor(n_estimators=100, max_depth=5)
        # regressor.fit(X_train, y_train)
        # y_pred = regressor.predict(X_test)
        rd = Ridge(alpha=alpha)
        rd.fit(X_train, y_train)
        y_pred = rd.predict(X_test)

        corr_list = []
        for ccc in range(0,132):
            correlation, _ = pearsonr(y_pred[:,ccc], y_test[:,ccc])
            corr_list.append(correlation)           
            # print("Pearson Correlation:", correlation)
        corrDF = np.hstack(corr_list)
        corrDF_DF = pd.DataFrame(corrDF.reshape(1,-1), columns=np.hstack([met_control.columns[2:], met_control.columns[2:]+"_stress"]))

        corrDF_DF.to_csv(f"{outputFolder}{modell}_cv{i}_alpha{alpha}.csv", index = False)
        # Pred_Corr.append(corrDF)
        

image_features_data = image_extractor(model_name=modell)
outputFolder = "../../temp/DLs_RR/"
for i in range(cv_1, cv_1+10):

    for alpha in [1, 100, 1000, 10000]:

        # log_file = open(outfile, 'a')
        # sys.stdout = log_file
        print(f"Now is running CV {i} RR alpha {alpha}")
        ll_reg(modell=modell,
                image_features_data = image_features_data, 
                i=i,
                alpha=alpha,
                outputFolder=outputFolder)
        # sys.stdout = sys.__stdout__  # Restore stdout
        # log_file.close()
    
