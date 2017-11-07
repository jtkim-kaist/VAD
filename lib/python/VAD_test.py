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
        opts, args = getopt.getopt(sys.argv[1:], 'hm:b:l:', ["data_dir=", "norm_dir=", "model_dir="])
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
        elif opt == '-b':
            batch_size = int(arg)
        elif opt == '-l':
            data_len = int(arg)
        elif opt == '--data_dir':
            data_dir = str(arg)
        elif opt == '--norm_dir':
            norm_dir = str(arg)
        elif opt == '--model_dir':
            model_dir = str(arg)

    if mode == 0:
        Vp.test_config(c_test_dir=data_dir,
                       c_norm_dir=norm_dir,
                       c_initial_logs_dir=model_dir, c_batch_size_eval=batch_size,
                       c_data_len=data_len)

        pred, label = Vp.main()

    elif mode == 1:
        Vb.test_config(c_test_dir=data_dir,
                       c_norm_dir=norm_dir,
                       c_initial_logs_dir=model_dir, c_batch_size_eval=batch_size,
                       c_data_len=data_len)

        pred, label = Vb.main()

    elif mode == 2:
        Vd.test_config(c_test_dir=data_dir,
                       c_norm_dir=norm_dir,
                       c_initial_logs_dir=model_dir, c_batch_size_eval=batch_size,
                       c_data_len=data_len)

        pred, label = Vd.main()

    elif mode == 3:

        Vl.test_config(c_test_dir=data_dir,
                       c_norm_dir=norm_dir,
                       c_initial_logs_dir=model_dir, c_batch_size_eval=batch_size,
                       c_data_len=data_len)

        pred, label = Vl.main()
    
    sio.savemat('./result/pred.mat', {'pred': pred})
    sio.savemat('./result/label.mat', {'label': label})
    print("done")
