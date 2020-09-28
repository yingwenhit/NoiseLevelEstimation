#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Power by Zongsheng Yue 2019-01-07 14:36:55
# for real image test

import os
import numpy as np
import scipy.io as scio
import sys
import cv2
from cv2 import imread
from utils import im2patch, im2double
import time

def noise_estimate(im, pch_size=8):
    '''
    Implement of noise level estimation of the following paper:
    Chen G , Zhu F , Heng P A . An Efficient Statistical Method for Image Noise Level Estimation[C]// 2015 IEEE International Conference
    on Computer Vision (ICCV). IEEE Computer Society, 2015.
    Input:
        im: the noise image, H x W x 3 or H x W numpy tensor, range [0,1]
        pch_size: patch_size
    Output:
        noise_level: the estimated noise level
    '''

    if im.ndim == 3:
        im = im.transpose((2, 0, 1))
    else:
        im = np.expand_dims(im, axis=0)
        # print('-------')

    # image to patch
    pch = im2patch(im, pch_size, 3)  # C x pch_size x pch_size x num_pch tensor
    num_pch = pch.shape[3]
    pch = pch.reshape((-1, num_pch))  # d x num_pch matrix
    d = pch.shape[0]

    mu = pch.mean(axis=1, keepdims=True)  # d x 1
    X = pch - mu
    sigma_X = np.matmul(X, X.transpose()) / num_pch
    sig_value, _ = np.linalg.eigh(sigma_X)
    sig_value.sort()

    for ii in range(-1, -d-1, -1):
        tau = np.mean(sig_value[:ii])
        if np.sum(sig_value[:ii]>tau) == np.sum(sig_value[:ii] < tau):
            return np.sqrt(tau)


if __name__ == '__main__':
    
    path = '/Users/wenying/OneDrive/0.PhD.Candidate/12. NPSID_IEEE/NPSID_IEEE/Matlab_code/RAISE-TestImage/'
    files = os.listdir(path)

    for file in files:
        if ".mat" in file:
            dataFile = os.path.join(path, file)
            data = scio.loadmat(dataFile)
            ImgName = ['I1', 'I2', 'I3']
            for Img in ImgName:
                im_noise = data[Img]
                
                start = time.time()
                est_level = noise_estimate(im_noise, 8)
                end = time.time()
                time_elapsed = end -start

                str_p = "Name: {0:s} Channel: {1:s} Time: {2:.4f}, Estimated Level: {3:6.4f}"
                print(str_p.format(file, Img, time_elapsed, est_level))



    # dataFile = 'Imat.mat'
    # data = scio.loadmat(dataFile)
    # im = data['I']
    # #im = imread(sys.argv[1],cv2.IMREAD_GRAYSCALE)
    # im = im2double(im)
    
    # im_noise = im

    # start = time.time()
    # est_level = noise_estimate(im_noise, 8)
    # end = time.time()
    # time_elapsed = end -start

    # str_p = "Time: {0:.4f}, Estimated Level: {1:6.4f}"
    # print(str_p.format(time_elapsed, est_level*255))


