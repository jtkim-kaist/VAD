import sys

sys.path.insert(0, './lib/python')
import VAD_Proposed as Vp
import VAD_DNN as Vd
import VAD_bDNN as Vb
import VAD_LSTM_2 as Vl
import scipy.io as sio
import os, getopt
# norm_dir = "./norm_data"
# data_dir = "./sample_data"
# ckpt_name = '/model9918and41.ckpt-2'
# model_dir = "./saved_model"
# valid_batch_size = 4134

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'h', ["data_dir=", "save_dir="])
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(1)

    if len(opts) != 2:
        print("arguments are not enough.")
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            sys.exit(0)
        elif opt == '--data_dir':
            data_dir = str(arg)
        elif opt == '--save_dir':
            save_dir = str(arg)

    data_dir = os.path.abspath('../..') + '/data' + data_dir
    train_data_dir = data_dir + '/train'
    valid_data_dir = data_dir + '/valid'

    save_dir = os.path.abspath('../..') + '/data' + save_dir
    train_save_dir = save_dir + '/train'
    valid_save_dir = save_dir + '/valid'

    os.system("rm -rf " + save_dir)
    os.system("mkdir " + save_dir)
    os.system("mkdir " + save_dir + '/train')
    os.system("mkdir " + save_dir + '/valid')

    os.system("matlab -r \"try acoustic_feat_ex(\'%s\',\'%s\'); catch; end; quit\"" % (train_data_dir, train_save_dir))
    os.system("matlab -r \"try acoustic_feat_ex(\'%s\',\'%s\'); catch; end; quit\"" % (valid_data_dir, valid_save_dir))

    # os.system("rm -rf")
    print("done")
