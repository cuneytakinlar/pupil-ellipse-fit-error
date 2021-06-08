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
    totalError = 0.0

    width = input_bk[0].shape[1]
    height = input_bk[0].shape[0]

    ellipseFitErrorDll.EllipseFitError.restype = ctypes.c_double
    count = 0
    for i in range(batchSize):        
        img = np.uint8(input_bk[i])
        img[img == 1] = 255

        gt = np.uint8(target_bk[i])
        gt[gt == 1] = 255
    
        error = ellipseFitErrorDll.EllipseFitError(width, height, ctypes.c_void_p(img.ctypes.data), ctypes.c_void_p(gt.ctypes.data))

        totalError += error
        if error != 0: count += 1

    if count == 0: return 0.0
    else: return totalError/count

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
    totalError = 0.0

    width = input_bk[0].shape[1]
    height = input_bk[0].shape[0]

    ellipseFitErrorDll.EllipseFitError.restype = ctypes.c_double
    count = 0
    for i in range(batchSize):        
        img = np.uint8(input_bk[i])
        img[img == 1] = 255
        
        gt = np.uint8(target_bk[i])
        gt[gt == 1] = 255

        error = ellipseFitErrorDll.EllipseFitError(width, height, ctypes.c_void_p(img.ctypes.data), ctypes.c_void_p(gt.ctypes.data))
        
        totalError += error
        if error != 0: count += 1

    if count == 0: return 0.0
    else: return totalError/count    
