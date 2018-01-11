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
        opts, args = getopt.getopt(sys.argv[1:], 'hm:e:', ["train_step=", "prj_dir="])
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(1)

    if len(opts) != 4:
        print("arguments are not enough.")
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            sys.exit(0)
        elif opt == '-m':
            mode = int(arg)
        elif opt == '-e':
            extract_feat = int(arg)
        elif opt == '--train_step':
            train_step = int(arg)
        elif opt == '--prj_dir':
            prj_dir = str(arg)

    data_dir = prj_dir + '/data/raw'
    print(data_dir)
    train_data_dir = data_dir + '/train'
    valid_data_dir = data_dir + '/valid'

    save_dir = prj_dir + '/data/feat'
    train_save_dir = save_dir + '/train'
    valid_save_dir = save_dir + '/valid'

    if extract_feat:

        os.system("rm -rf " + save_dir)
        os.system("mkdir " + save_dir)
        os.system("mkdir " + save_dir + '/train')
        os.system("mkdir " + save_dir + '/valid')
        os.system(
            "matlab -r \"try acoustic_feat_ex(\'%s\',\'%s\'); catch; end; quit\"" % (train_data_dir, train_save_dir))
        os.system(
            "matlab -r \"try acoustic_feat_ex(\'%s\',\'%s\'); catch; end; quit\"" % (valid_data_dir, valid_save_dir))

        train_norm_dir = save_dir + '/train/global_normalize_factor.mat'
        test_norm_dir = prj_dir + '/norm_data/global_normalize_factor.mat'

        os.system("cp %s %s" % (train_norm_dir, test_norm_dir))

    if mode == 0:

        logs_dir = prj_dir + '/logs'
        os.system("rm -rf " + logs_dir + '/train')
        os.system("rm -rf " + logs_dir + '/valid')
        os.system("mkdir " + logs_dir + '/train')
        os.system("mkdir " + logs_dir + '/valid')

        Vd.train_config(save_dir+'/train', save_dir+'/valid', prj_dir+'/logs', 256,
                        train_step, 'train')

        Vd.main()
    # os.system("rm -rf")
    print("done")