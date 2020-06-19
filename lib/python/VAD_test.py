import sys
sys.path.insert(0, './lib/python')
import VAD_Proposed as Vp
import VAD_DNN as Vd
import VAD_bDNN as Vb
import VAD_LSTM_2 as Vl
import scipy.io as sio
import graph_test as graph_test
import os, getopt
import glob

from time import time

# norm_dir = "./norm_data"
# data_dir = "./sample_data"
# ckpt_name = '/model9918and41.ckpt-2'
# model_dir = "./saved_model"
# valid_batch_size = 4134

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

    if mode == 0:
        # Vp.test_config(c_test_dir=data_dir,
        #                c_norm_dir=norm_dir,
        #                c_initial_logs_dir=model_dir, c_batch_size_eval=batch_size,
        #                c_data_len=data_len)
        #
        # pred, label = Vp.main()

        if is_default:
            graph_list = sorted(glob.glob(model_dir + '/backup/backup_pb/frozen_model_ACAM.pb'))
            norm_dir = model_dir + '/backup/backup_norm'
            pred, label = graph_test.do_test(graph_list[-1], data_dir, norm_dir, data_len, is_default, mode)
        else:
            graph_list = sorted(glob.glob(model_dir + '/graph/ACAM/*.pb'))
            print(graph_list)
            pred, label = graph_test.do_test(graph_list[-1], data_dir, norm_dir, data_len, is_default, mode)

    elif mode == 1:

        print(os.path.abspath('./configure/bDNN'))
        sys.path.insert(0, os.path.abspath('./configure/bDNN'))

        import config as cg
        if is_default:
            graph_list = sorted(glob.glob(model_dir + '/backup/backup_pb/frozen_model_bDNN.pb'))
            norm_dir = model_dir + '/backup/backup_norm'
            pred, label = graph_test.do_test(graph_list[-1], data_dir, norm_dir, data_len, is_default, mode)
        else:
            graph_list = sorted(glob.glob(model_dir + '/graph/bDNN/*.pb'))
            print(graph_list)
            pred, label = graph_test.do_test(graph_list[-1], data_dir, norm_dir, data_len, is_default, mode)

        # Vb.test_config(c_test_dir=data_dir,
        #                c_norm_dir=norm_dir,
        #                c_initial_logs_dir=model_dir, c_batch_size_eval=batch_size,
        #                c_data_len=data_len)
        # Vb.test_config(c_test_dir=data_dir,
        #                c_norm_dir='/home/sbie/storage3/github/VAD_Toolkit/VAD/saved_model/backup_norm',
        #                c_initial_logs_dir='/home/sbie/storage3/github/VAD_Toolkit/VAD/saved_model/backup_ckpt', c_batch_size_eval=batch_size,
        #                c_data_len=data_len)
        #
        # pred, label = Vb.main()

    elif mode == 2:

        start_time = time()
        print(os.path.abspath('./configure/DNN'))
        sys.path.insert(0, os.path.abspath('./configure/DNN'))

        import config as cg

        if is_default:
            graph_list = sorted(glob.glob(model_dir + '/backup/backup_pb/frozen_model_DNN.pb'))
            norm_dir = model_dir + '/backup/backup_norm'
            pred, label = graph_test.do_test(graph_list[-1], data_dir, norm_dir, data_len, is_default, mode)
        else:
            graph_list = sorted(glob.glob(model_dir + '/graph/DNN/*.pb'))
            pred, label = graph_test.do_test(graph_list[-1], data_dir, norm_dir, data_len, is_default, mode)

        # Vd.test_config(c_test_dir=data_dir,
        #                c_norm_dir=norm_dir,
        #                c_initial_logs_dir=model_dir, c_batch_size_eval=batch_size,
        #                c_data_len=data_len)
        #
        # pred, label = Vd.main()

        end_time = time()

        time_taken = end_time - start_time
        print(time_taken)

    elif mode == 3:

        sys.path.insert(0, os.path.abspath('./configure/LSTM'))

        import config as cg

        if is_default:
            graph_list = sorted(glob.glob(model_dir + '/backup/backup_pb/frozen_model_LSTM.pb'))
            norm_dir = model_dir + '/backup/backup_norm'
            pred, label = graph_test.do_test(graph_list[-1], data_dir, norm_dir, data_len, is_default, mode)
        else:
            graph_list = sorted(glob.glob(model_dir + '/graph/LSTM/*.pb'))
            print(graph_list)
            pred, label = graph_test.do_test(graph_list[-1], data_dir, norm_dir, data_len, is_default, mode)

        # Vl.test_config(c_test_dir=data_dir,
        #                c_norm_dir=norm_dir,
        #                c_initial_logs_dir=model_dir, c_batch_num=200, c_seq_size=20,
        #                c_data_len=data_len)
        #
        # pred, label = Vl.main()
    
    sio.savemat('result/pred.mat', {'pred': pred})
    sio.savemat('result/label.mat', {'label': label})
    print("done")
