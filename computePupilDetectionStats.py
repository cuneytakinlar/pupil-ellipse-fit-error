import os
import numpy as np
import math
from os import listdir
from os.path import isfile, join

users = [
    [],         # dir = 0 (non-existent)
    [1, 4, 9],  # dir = 1
    [4, 10, 13],  # dir = 2
    [16, 19, 21],  # dir = 3
    [1, 2, 12],  # dir = 4
    [6, 10, 11],  # dir = 5
    [2, 5, 13],  # dir = 6
    [15, 18, 21],  # dir = 7
    [2, 7, 9],  # dir = 8
    [16, 17, 18],  # dir = 9
    [1, 8, 11],  # dir = 10
    [2, 7, 13],  # dir = 11
    [1, 2, 9],  # dir = 12
    [1, 2, 9],  # dir = 13
    [10, 17, 22],  # dir = 14
    [1, 2, 7],  # dir = 15
    [1, 2, 13],  # dir = 16
    [3, 5, 12],  # dir = 17
    [2, 7, 11],  # dir = 18
    [2, 3, 6],  # dir = 19
    [3, 4, 7],  # dir = 20
    [4, 11, 12],  # dir = 21
    [1, 2, 17],  # dir = 22
]

gt_path = './LPW'
result_path = './result_files'

def ComputeStatForOneFile(dir, maxDistance):
    result = np.zeros(maxDistance)
    noPupils = 0

    for user in users[dir]:
        # Read the gt
        gtFile = f"{gt_path}/{dir}/{user}.txt"
        contents = open(gtFile, "r").read().split("\n")
        gt = []
        for line in contents:
            if line == "":
                break

            x, y = line.split(' ')
            gt.append((float(x), float(y)))

        # Read the predictions
        predFile = f"{result_path}/{dir}/{user}.txt"
        contents = open(predFile, "r").read().split("\n")
        pred = []
        for line in contents:
            if line == "":
                break

            x, y = line.split(' ')
            pred.append((float(x), float(y)))

#        print(len(pred))
        if len(pred) < len(gt):
            print(f"Length of pred {len(pred)} < the length of gt {len(gt)}!")
            os._exit(0)

        # Compute the detection accuracies as the distance is increased
        noPupils += len(pred)
        for index in range(len(pred)):
            x1, y1 = gt[index]
            x2, y2 = pred[index]

            dx = x1 - x2
            dy = y1 - y2
            distance = math.sqrt(dx*dx + dy*dy)
            for d in range(maxDistance):
                if distance <= d:
                    result[d] += 1

    return noPupils, result

#---------------------------------------------------------------------
# Main
#
if __name__ == '__main__':    
    maxDistance = 15+1

    results = []
    totalNoImgs = 0

    for dir in range(1, 23, 1):
        noImgs, result = ComputeStatForOneFile(dir, maxDistance)
        results.append((noImgs, result));        
        totalNoImgs += noImgs
        print(f"{100*result[5]/noImgs:6.2f}")

    sum = np.zeros(maxDistance)
    for index in range(maxDistance):
        for j in range(len(results)):
            sum[index] += results[j][1][index]                   # Weighted average 

    for index in range(maxDistance):
        print(f"{index:3d} {100*sum[index]/totalNoImgs:6.2f}")   # Weighted average 
