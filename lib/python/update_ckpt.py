import sys
sys.path.insert(0, './lib/python')
import VAD_Proposed as Vp
import VAD_DNN as Vd
import VAD_bDNN as Vb
import VAD_LSTM_2 as Vl
import scipy.io as sio
import os, getopt
import glob

# norm_dir = "./norm_data"
# data_dir = "./sample_data"
# ckpt_name = '/model9918and41.ckpt-2'
# model_dir = "./saved_model"
# valid_batch_size = 4134

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hu:', ["model=", "prj_dir="])

    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(1)

    if len(opts) != 3:
        print("arguments are not enough.")
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            sys.exit(0)
        if opt == '-u':
            update = int(arg)
        elif opt == '--model':
            model = str(arg)
        elif opt == '--prj_dir':
            prj_dir = str(arg)

    logs_dir = prj_dir + '/logs'

    if update == 1:
        ckpt_list = sorted(glob.glob(logs_dir + '/model*'))
        new_ckpt = old_ckpt= ckpt_list[-3:]
        new_ckpt = [(ckpt_name.replace('model', model)) for ckpt_name in new_ckpt]
        save_dir = [prj_dir + '/saved_model/' + os.path.basename(ckpt_name) for ckpt_name in new_ckpt]
        save_dir = [ckpt_name.split('.')[0] + '.' + ckpt_name.split('.')[-1] for ckpt_name in save_dir]

        for x, y, z in zip(old_ckpt, new_ckpt, save_dir):
            os.system('cp -f ' + x + ' ' + y)
            os.system('mv -f ' + y + ' ' + z)
        print("checkpoint update done!")
    else:
        backup_dir = prj_dir + '/saved_model/backup_ckpt'
        backup_ckpt = sorted(glob.glob(backup_dir + '/' + model + '*'))
        restore_ckpt = [prj_dir + '/saved_model/' + os.path.basename(ckpt_name) for ckpt_name in backup_ckpt]

        for x, y in zip(backup_ckpt, restore_ckpt):
            os.system('cp -f ' + x + ' ' + y)

        backup_norm_dir = prj_dir + '/norm_data/backup_norm/global_normalize_factor.mat'
        norm_dir = prj_dir + '/norm_data/global_normalize_factor.mat'

        os.system('cp -f ' + backup_norm_dir + ' ' + norm_dir)

        print("checkpoint restore done!")
