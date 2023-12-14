from tensorflow.keras.applications.vgg16 import VGG16
from tensorflow.keras.applications.vgg16 import preprocess_input, decode_predictions

import numpy as np

model = VGG16(weights='imagenet')


img = np.random.randint(0, 256, size=(224, 224, 3), dtype=np.uint8)
# img = img.img_to_array(img)
img = np.expand_dims(img, axis=0)
img = preprocess_input(img)
temp = model.predict(img)
print(temp)
