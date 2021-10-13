import os
import cv2
import tensorflow as tf
import sys
import random

Size_X, Size_Y = 320, 240
#Size_X, Size_Y = 640, 480

def _bytes_feature(value):
    return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))

# Directory where the datasets will be written
output_path = './tfRecords'
if not os.path.exists(output_path): os.mkdir(output_path)

# Create TF writers
trainWriter = tf.io.TFRecordWriter(f"{output_path}/{Size_X}x{Size_Y}-train")
validationWriter = tf.io.TFRecordWriter(f"{output_path}/{Size_X}x{Size_Y}-validation")

# LPW GT data
grayImgFolderPath = "./LPW/gt/imgs/"
GTImgFolderPath = "./LPW/gt/masks/"

# Get all files & shuffle them
filenamelist = os.listdir(GTImgFolderPath)
random.shuffle(filenamelist) 

# Split ratio
splitRatio = 0.8    # 80% training, 20% validation

# Number of training images
noOfTrainingImages = int(len(filenamelist)*splitRatio)
noOfValidationImages = len(filenamelist) - noOfTrainingImages

# Number of augmentations
noOfAugmentations = 5   # Can be increased if desired

image_no = 0
for num_split_dataset in range(noOfAugmentations):    
    for j in range(len(filenamelist)):  # For each file
        image_no = image_no + 1
        if image_no % 100 == 0: print(f"{image_no} images completed...")

        grayimg_name = grayImgFolderPath + '/' + filenamelist[j]
        GTimg_name = GTImgFolderPath + '/' + filenamelist[j]

        gray_img = cv2.imread(grayimg_name, cv2.IMREAD_GRAYSCALE)
        binaryGT_img = cv2.imread(GTimg_name, cv2.IMREAD_GRAYSCALE)
        height, width = binaryGT_img.shape
        random_degree = random.uniform(-5.0, 5.0)
        random_scale = random.uniform(0.95, 1.05)
        matrix = cv2.getRotationMatrix2D((width/2, height/2), random_degree, random_scale)

        ########################################## make Affine
        gray_img = cv2.warpAffine(gray_img, matrix, (width, height))
        binaryGT_img = cv2.warpAffine(binaryGT_img, matrix, (width, height))
        binaryGT_img = cv2.threshold(binaryGT_img, 127, 255, cv2.THRESH_BINARY)[1]
        ##################################################

        if gray_img.shape[1] != Size_X:
            gray_img = cv2.resize(gray_img, (Size_X, Size_Y), interpolation=cv2.INTER_CUBIC)
            binaryGT_img = cv2.resize(binaryGT_img, (Size_X, Size_Y), interpolation=cv2.INTER_CUBIC)
            binaryGT_img = cv2.threshold(binaryGT_img, 127, 255, cv2.THRESH_BINARY)[1]

        # Create a feature
        feature = {'train/image': _bytes_feature(tf.compat.as_bytes(gray_img.tobytes())),
                    'train/label': _bytes_feature(tf.compat.as_bytes(binaryGT_img.tobytes()))
                  }
        # Create an example protocol buffer
        example = tf.train.Example(features=tf.train.Features(feature=feature))

        if j < noOfTrainingImages: trainWriter.write(example.SerializeToString())
        else: validationWriter.write(example.SerializeToString())

print(f"Number of training images in the dataset: {noOfTrainingImages*noOfAugmentations}")
print(f"Number of validation images in the dataset: {noOfValidationImages*noOfAugmentations}")
print(f"Total number of images: {image_no}")

trainWriter.close()
validationWriter.close()




