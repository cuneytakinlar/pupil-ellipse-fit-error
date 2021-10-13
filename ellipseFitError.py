import torch
import numpy as np

from ctypes import *
import ctypes

ellipseFitErrorDll = CDLL("./EllipseFitError.dll")

#--------------------------------------------------------------------------
# This is used with UNet
#
def EllipseFitErrorUNet(input, target):
    input_bk = input[:, 0].clone().detach().cpu().numpy()
    input_bk[input_bk < 0.5] = 0
    input_bk[input_bk >= 0.5] = 1

    target_bk = target[:, 0].clone().detach().cpu().numpy()
    target_bk[target_bk < 0.5] = 0
    target_bk[target_bk >= 0.5] = 1

    batchSize = input_bk.shape[0]

    width = input_bk[0].shape[1]
    height = input_bk[0].shape[0]

    ellipseFitErrorDll.EllipseFitError.restype = ctypes.c_double
    distances = np.zeros(width*height, np.float64)
    totalNoPixels = 0
    sum = torch.tensor(0.0) 
    for i in range(batchSize):        
        img = np.uint8(input_bk[i])
        img[img == 1] = 255

        gt = np.uint8(target_bk[i])
        gt[gt == 1] = 255
    
        noPixels = ellipseFitErrorDll.EllipseFitError(width, height, ctypes.c_void_p(img.ctypes.data), ctypes.c_void_p(gt.ctypes.data), ctypes.c_void_p(distances.ctypes.data))
        if noPixels == 0: continue

        sum += torch.sum(torch.from_numpy(np.resize(distances, int(noPixels))))        
        totalNoPixels += noPixels

    if totalNoPixels == 0: return sum
    else:
        sum = torch.tensor(sum.item(), requires_grad=True)
        return torch.sqrt(sum/totalNoPixels)      

#--------------------------------------------------------------------------
# This is used with DenseNet
#
def EllipseFitErrorDenseNet(input, target):
    input_bk = input.clone().detach().cpu().numpy()
    input_bk[input_bk < 0.5] = 0
    input_bk[input_bk >= 0.5] = 1

    target_bk = target.clone().detach().cpu().numpy()
    target_bk[target_bk < 0.5] = 0
    target_bk[target_bk >= 0.5] = 1

    batchSize = input_bk.shape[0]

    width = input_bk[0].shape[1]
    height = input_bk[0].shape[0]

    ellipseFitErrorDll.EllipseFitError.restype = ctypes.c_double
    distances = np.zeros(width*height, np.float64)
    totalNoPixels = 0
    sum = torch.tensor(0.0)
    for i in range(batchSize):        
        img = np.uint8(input_bk[i])
        img[img == 1] = 255
        
        gt = np.uint8(target_bk[i])
        gt[gt == 1] = 255

        noPixels = ellipseFitErrorDll.EllipseFitError(width, height, ctypes.c_void_p(img.ctypes.data), ctypes.c_void_p(gt.ctypes.data), ctypes.c_void_p(distances.ctypes.data))
        if noPixels == 0: continue

        sum += torch.sum(torch.from_numpy(np.resize(distances, int(noPixels))))        
        totalNoPixels += noPixels

    if totalNoPixels == 0: return sum
    else:
        sum = torch.tensor(sum.item(), requires_grad=True) 
        return torch.sqrt(sum/totalNoPixels)