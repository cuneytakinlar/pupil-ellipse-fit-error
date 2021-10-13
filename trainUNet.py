from __future__ import print_function

import torch
import os
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np
import tensorflow as tf
from torchvision.utils import save_image

import cv2
import warnings
import math

from unet import UNet, MyUNet
from ellipseFitError import EllipseFitErrorUNet

#----------------------- _parse_image_function ----------------------
feature = {'train/image': tf.io.FixedLenFeature([], tf.string),
           'train/label': tf.io.FixedLenFeature([], tf.string)}

# Parse the input tf.Example proto using the dictionary above.
def _parse_image_function(example_proto):
  return tf.io.parse_single_example(example_proto, feature)

#----------------------- processData ----------------------
# Given a model & dataset, process the dataset on the model using the processing mode
#
def processData(model, dataSet, epoch, processingMode = "train"):
    global imageNo
    imageNo = 1

    # Access the data in batches of batchSize
    dataSetByBatch = dataSet.batch(batch_size=batchSize)

    noImages = 0
    current_batch = 0
    total_loss = 0
    interval_loss = 0
    for image_features in dataSetByBatch:            
        image = tf.io.decode_raw(image_features['train/image'], tf.uint8)
        label = tf.io.decode_raw(image_features['train/label'], tf.uint8)

        image = tf.reshape(image, [-1, Size_Y, Size_X])
        image = tf.expand_dims(image, 1)

        label = tf.reshape(label, [-1, Size_Y, Size_X])
        label = tf.expand_dims(label, 1)

        noImages += len(image)
        current_batch += 1

        # Convert the images to [0-1] & load them to GPU
        image = image.numpy().astype(np.float32)/255
        label = label.numpy().astype(np.float32)/255

        image = torch.from_numpy(image)
        image = image.to(device)

        label = torch.from_numpy(label)
        label = label.to(device)

        # Pass the image through the model & then sigmoid
        output = model(image)
        output = nn.functional.sigmoid(output)

        # Calculate loss
        loss = F.binary_cross_entropy(output, label)

        # Add EllipseFitError regularization Term
        ellipseFitError = EllipseFitErrorUNet(output, label)*0.1
        print(ellipseFitError.item())
        loss += ellipseFitError

        if processingMode == "train":
            # Do backward propagation in training mode & update weights
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        total_loss += float(loss)
        interval_loss += float(loss)

        # Print information on the screen
        DIV = 100
        if current_batch % DIV == 0:
            avg_loss = total_loss/current_batch
            avg_interval_loss = interval_loss/DIV
            interval_loss = 0
            print(f"Epoch: {epoch}, Batch: {current_batch:4}, avg. loss: {avg_loss:.5f}, int. loss: {avg_interval_loss:.5f}")

            if processingMode == "train":
                fileno = int(current_batch/DIV)
                filename = f"output/{fileno:05d}.png"
                save_image(output[0], filename)               

        del image, label, output, loss

    return total_loss/current_batch # Return the average loss during the epoch 

#-------------------------- main function -------------------------------
def main():
    # Make the following variables global so that other function can access them
    global Size_X, Size_Y, batchSize, device, optimizer

    torch.manual_seed(1)
    warnings.filterwarnings("ignore")

    use_cuda = torch.cuda.is_available()
    device = torch.device("cuda" if use_cuda else "cpu")
    print("DEVICE: ", device)

    shuffle_buffer = int(2000)
    Size_X, Size_Y, batchSize = 320, 240, 16
#    Size_X, Size_Y, batchSize = 640, 480, 10
    
    # Datasets
    rawTrainingDataset = tf.data.TFRecordDataset(f"./tfRecords/{Size_X}x{Size_Y}-train")
    rawValidationDataset = tf.data.TFRecordDataset(f"./tfRecords/{Size_X}x{Size_Y}-validation")

    trainingSet = rawTrainingDataset.map(_parse_image_function)
    validationSet = rawValidationDataset.map(_parse_image_function)

    # This is where we will save the trained models at the end of each epoch
    trained_model_path = './trained_model'
    if not os.path.exists(trained_model_path): os.mkdir('./trained_model')

    # This is where we will output debug output images
    output_path = './output'
    if not os.path.exists(output_path): os.mkdir(output_path)

    # Create a UNet model
    learning_rate = 1e-4
    model = MyUNet(32)

    # Load a previous saved model (if it exists)
#    model_file_name = "trained_model/model1.pt"    
#    model.load_state_dict(torch.load(model_file_name))

    model.to(device)    # Upload the model to GPU    

    # Create an optimizer
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', patience=2)

    # Open log file
    logFile = open("log.txt", "w")

    noEpochs = 30    
    for epoch in range(1, noEpochs+1, 1):
        print(f"==============Training Epoch: {epoch}===============")

        #----------------- Training ------------------------
        # Shuffle training data & get in batches
        # Put the model in train mode
        model.train()   
        torch.set_grad_enabled(True) # Explicitly enable gradient computation    

        trainingSet = trainingSet.shuffle(buffer_size=shuffle_buffer)
        avg_training_loss = processData(model, trainingSet, epoch, "train")
        print(f"----> End of training epoch {epoch:2d}. avg. loss: {avg_training_loss:.5f}")

        #----------------- Validation ------------------------
        # Put the model in evaluation mode & do validation
        model.eval() 
        torch.set_grad_enabled(False) # Explicitly disable gradient computation

        avg_val_loss = processData(model, validationSet, epoch, "validation")
        print(f"----> End of validation. epoch {epoch:2d}. avg. loss: {avg_val_loss:.5f}")
        scheduler.step(avg_val_loss)

        # Write the average loss for the epoch to the log file
        logFile.write(f"{epoch:4d}{avg_training_loss:10.5f}{avg_val_loss:10.5f}\n")
        logFile.flush()

    # Save the trained model to disk
    savename = f"{trained_model_path}/Unet-trained-model.pt"
    torch.save(model.state_dict(), savename)

    logFile.close()

#-------------------------------------------------------
# Call the main function & start processing
if __name__ == '__main__':
    main()