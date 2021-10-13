from __future__ import print_function

import torch
import os
import sys
import numpy as np
import cv2
from unet import MyUNet

def main():
#    Size_X, Size_Y = 320, 240
    Size_X, Size_Y = 640, 480

    MULT = 640/Size_X

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("DEVICE: ", device)

    model = MyUNet(32)
#    model_file_name = "trained_model/modelXX.pt"    
    model_file_name = "efe-Unet-trained-model-640x480.pt" # This is the model for ellipse fit error
    model.load_state_dict(torch.load(model_file_name))

    model.eval()
    model.to(device)
    print('Test Started!')

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

    image_path = './LPW'
    result_path = './result_files'
    if not os.path.exists(result_path):
        os.mkdir(result_path)

    no_images_with_more_than_one_blob = 0

    for dir in range(1, 23, 1):
        for user in users[dir]:
            resultDir = f"{result_path}/{dir}"
            if not os.path.exists(resultDir):
                os.mkdir(resultDir)

            resultTxt = open(f"{result_path}/{dir}/{user}.txt", 'w')

            aviFile = f"{image_path}/{dir}/{user}.avi"
            print(aviFile)
            
            prev_x, prev_y = 0, 0
            frameNo = 0
            cap = cv2.VideoCapture(aviFile)
            while True:
                ret, inputImg = cap.read()     
                if ret == False: break           

                inputImg = inputImg[:, :, 0]  # Take the first plane
                inputImg = cv2.resize(inputImg, (Size_X, Size_Y), interpolation=cv2.INTER_CUBIC)
                inputImg = inputImg[np.newaxis, np.newaxis, :]
                inputImg = inputImg.astype(np.float32)/255     

                image = torch.from_numpy(inputImg)
                image = image.to(device)
                output = model(image)
                output_bk = output[:, 0].clone().detach().cpu().numpy()
                ttt = output_bk
                ttt[ttt < 0.5] = 0
                ttt[ttt >= 0.5] = 1
                if np.count_nonzero(ttt) == 0:
                    output_bk[output_bk < 0.25] = 0
                    output_bk[output_bk >= 0.25] = 1
                else:
                    output_bk[output_bk < 0.5] = 0
                    output_bk[output_bk >= 0.5] = 1

                ## Connected Component Analysis
                if np.count_nonzero(output_bk) != 0:
                    no_labels, labels, stats, center = cv2.connectedComponentsWithStats(output_bk[0, :, :].astype(np.uint8))

                    if no_labels > 2:
                        no_images_with_more_than_one_blob += 1

                    stats = stats[1:, :]
                    pupil_candidate = np.argmax(stats[:, 4]) + 1
                    centerX = center[pupil_candidate][0]*MULT
                    centerY = center[pupil_candidate][1]*MULT

                    prev_x = centerX
                    prev_y = centerY
                    txt = f"{round(centerX, 3):.3f} {round(centerY, 3):.3f}\n"

                else:
                    txt = f"{round(prev_x, 3):.3f} {round(prev_y, 3):.3f}\n"
                resultTxt.write(txt)

    #            savename = result_path + file_name
    #            cv2.imwrite(savename, (output_bk[0]*255+inputImg_BK)/2)

                frameNo += 1
                if frameNo % 500 == 0:
                    print(f"frameNo: {frameNo}, no_images_with_more_than_one_blob: {no_images_with_more_than_one_blob}");

            resultTxt.close()
            sys.stdout.flush()

if __name__ == '__main__':
    main()
