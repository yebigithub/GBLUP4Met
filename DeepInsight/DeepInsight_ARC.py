##failed because cannot use the function "robjects.r" in ARC.

from pyDeepInsight import ImageTransformer
from pyDeepInsight.utils import Norm2Scaler
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.manifold import TSNE
# from sklearn.decomposition import NMF
# import umap.umap_ as umap

import pandas as pd
import numpy as np
import pyreadr

# Load the Rdata file

def SNP_img_generator(chrr="1", resolution=277):

    # Get the list from the Rdata file
    result = pyreadr.read_r('../../Geno/SNPchrs/SNPchr'+ chrr+ '.RData')
    geno_chr = result["geno_chr"]

    expr = geno_chr

    y = expr.index
    X = expr.values
    X_train = X
    ln = Norm2Scaler()
    X_train_norm = ln.fit_transform(X_train)
    # X_test_norm = ln.transform(X_test)
    # Create t-SNE object
    distance_metric = 'cosine'
    reducer = TSNE(
        n_components=2,
        metric=distance_metric,
        init='random',
        learning_rate='auto',
        n_jobs=-1,
        random_state=42
    )

    # Initialize image transformer.
    pixel_size = (resolution,resolution)
    it = ImageTransformer(
        feature_extractor=reducer, 
        pixels=pixel_size)

    # Train image transformer on training data and transform training 
    # and testing sets. Values should be between 0 and 1.
    it.fit(X_train, plot=False)

    X_train_img = it.transform(X_train_norm, empty_value=1)
    # X_test_img = it.transform(X_test_norm, empty_value=1)
    np.save(f'../../SNPimg/Chr{Chr}_tsne_{resolution}.npy', X_train_img)

    print(f'Done with Chr {Chr}, Res {resolution}')



for Chr in range(2, 13):
        Chr = str(Chr)
        res = 277
        SNP_img_generator(chrr=Chr, resolution=res)



