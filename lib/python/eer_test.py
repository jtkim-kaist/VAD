import numpy as np
import scipy.io as sio
import sys
import os, sys, getopt
from sklearn import metrics
from scipy.optimize import brentq
from scipy.interpolate import interp1d


def eer(pred, label):

    fpr, tpr, thresholds = metrics.roc_curve(label, pred, pos_label=1)
    # valid_auc = metrics.auc(fpr, tpr)
    eer_result = brentq(lambda x: 1. - x - interp1d(fpr, tpr)(x), 0., 1.)
    return eer_result * 100


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'h', ["data_dir="])
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(1)

    # if len(opts) != 6:
    #     print("arguments are not enough.")
    #     sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            sys.exit(0)

        elif opt == '--data_dir':
            data_dir = str(arg)

    label = data_dir + '/label.mat'
    pred = data_dir + '/pred.mat'

    label = sio.loadmat(label)
    label = label['label']
    pred = sio.loadmat(pred)
    pred = pred['pred']

    eer_result = eer(pred, label)
    print(eer_result)



if __name__ == '__main__':
    main()