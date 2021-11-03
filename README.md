# Accurate CNN-based Pupil Segmentation with an Ellipse Fit Error Regularization Term

Cuneyt Akinlar, Hatice Kubra Kucukkartal, Cihan Topal

## Abstract

Semantic segmentation of images by Fully Convolutional Neural Networks (FCN) has gained increased attention in recent years as FCNs greatly outperform traditional segmentation algorithms. In this paper we propose using {\em Ellipse Fit Error} as a shape prior regularization term that can be added to a pixel-wise loss function, e.g., binary cross entropy, to train a CNN for pupil segmentation. We evaluate the performance of the proposed method by training a lightweight UNet architecture, and use three widely used real-world datasets for pupil center estimation, i.e., ExCuSe, ElSe, and Labeled Pupils in the Wild (LPW), containing a total of $\sim$230.000 images for performance evaluation. Experimental results show that the proposed method gives the best-known pupil detection rates for all datasets.

## Paper & Citation

This paper has been accepted for publication in Elsevier Expert Systems with Applications. If you use the code given here, please cite this paper <br>
Cuneyt Akinlar, Hatice Kubra Kucukkartal, Cihan Topal. Accurate CNN-based Pupil Segmentation with an Ellipse Fit Error Regularization Term. Expert Systems with Applications, vol. 188, pp. 116004, Feb. 2022. https://doi.org/10.1016/j.eswa.2021.116004.

## LPW Dataset & GT Segmentation Maps

Full LPW dataset can be downloaded from http://datasets.d2.mpi-inf.mpg.de/tonsen/LPW.zip <br>
GT segmentation maps for the LPW dataset can be downloaded from  https://www.kaggle.com/cuneytakinlar/lpw-gt-segmentation-maps <br>

## Installing the required Python packages

Create a virtual environment (we used conda). When you activate that virtual environment, type the following to download and install the required packages: <br>

% conda install tensorflow opencv torchvision numpy cudatoolkit -c pytorch <br>

## How to run the code

Before you run the code, you must select the image size to be used in dataset creation, training & test. Inside each file there is a line that you need to change. By default, these values are set to width=320, height=240. You can set them to any value that you want.

(1) Create the training & validatin datasets: <br>

% python createTFData.py <br>

(2) Train the model. There are two semantic segmentation models that you can use: UNet & DenseNet. You can train any one of them that you like. As the training continues, current model files are saved under trained_model directory after each epoch. <br>

% python trainUNet.py <br>

(3) Test the model. Edit testUNet.py & type in the name of the model file to use. Then: <br>

% python testUNet.py <br>

(4) Compute the pupil detection rates: <br>

% python computePupilDetectionStats.py <br>

Alternatively, you can use DenseNet for training & testing. In that case you need to use trainDenseNet.py & testDenseNet.py.










