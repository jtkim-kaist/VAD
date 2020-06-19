import sys
sys.path.insert(0, './lib/python')
import scipy.io as sio

from . import graph_test
import os, getopt
import glob

from time import time

# norm_dir = "./norm_data"
# data_dir = "./sample_data"
# ckpt_name = '/model9918and41.ckpt-2'
# model_dir = "./saved_model"
# valid_batch_size = 4134

def get_predictions(mode, data_len, is_default, data_dir, norm_dir, model_dir):
    if mode == 0:
        if is_default:
            graph_list = sorted(glob.glob(model_dir + '/backup/backup_pb/frozen_model_ACAM.pb'))
            norm_dir = model_dir + '/backup/backup_norm'
        else:
            graph_list = sorted(glob.glob(model_dir + '/graph/ACAM/*.pb'))
            print(graph_list)

    elif mode == 1:

        print(os.path.abspath('./configure/bDNN'))
        sys.path.insert(0, os.path.abspath('./configure/bDNN'))

        if is_default:
            graph_list = sorted(glob.glob(model_dir + '/backup/backup_pb/frozen_model_bDNN.pb'))
            norm_dir = model_dir + '/backup/backup_norm'
        else:
            graph_list = sorted(glob.glob(model_dir + '/graph/bDNN/*.pb'))
            print(graph_list)

    elif mode == 2:

        start_time = time()
        print(os.path.abspath('./configure/DNN'))
        sys.path.insert(0, os.path.abspath('./configure/DNN'))

        if is_default:
            graph_list = sorted(glob.glob(model_dir + '/backup/backup_pb/frozen_model_DNN.pb'))
            norm_dir = model_dir + '/backup/backup_norm'
        else:
            graph_list = sorted(glob.glob(model_dir + '/graph/DNN/*.pb'))

        end_time = time()

        time_taken = end_time - start_time
        print(time_taken)

    elif mode == 3:
        sys.path.insert(0, os.path.abspath('./configure/LSTM'))
        if is_default:
            graph_list = sorted(glob.glob(model_dir + '/backup/backup_pb/frozen_model_LSTM.pb'))
            norm_dir = model_dir + '/backup/backup_norm'
        else:
            graph_list = sorted(glob.glob(model_dir + '/graph/LSTM/*.pb'))
            print(graph_list)
    pred, label = graph_test.do_test(graph_list[-1], data_dir, norm_dir, data_len, is_default, mode)
    print('{} pred : {} label'.format(pred, label))
    return pred, label


if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hm:l:d:', ["data_dir=", "norm_dir=", "model_dir="])
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(1)

    if len(opts) != 6:
        print("arguments are not enough.")
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            sys.exit(0)
        elif opt == '-m':
            mode = int(arg)
        elif opt == '-l':
            data_len = int(arg)
        elif opt == '-d':
            is_default = int(arg)
        elif opt == '--data_dir':
            data_dir = str(arg)
        elif opt == '--norm_dir':
            norm_dir = str(arg)
        elif opt == '--model_dir':
            model_dir = str(arg)

    pred, label = get_predictions(mode, data_len, is_default, data_dir, norm_dir, model_dir)
    
    sio.savemat('result/pred.mat', {'pred': pred})
    sio.savemat('result/label.mat', {'label': label})
    print("done")
