import os

import librosa
import numpy as np
import scipy.io as sio

from . import mrcg as mrcg

def th_classifier(input, th):
    output = (input >= th).choose(input, 1)
    output = (output < th).choose(output, 0)
    return output


def binary_saver(name_file, data, num_file ):
    bin_num = '%3.3d'%num_file
    mrcg_name = name_file + '_'+bin_num+'.bin'
    spec_name = name_file + '_spec_' + bin_num+'.txt'
    fid_file = open(mrcg_name,'wb')
    fid_txt = open(spec_name,'wt')
    data.tofile(fid_file)
    nc, nr = np.shape(data)
    fid_txt.write(str(nc)+','+str(nr)+',float32')
    fid_file.close()
    fid_txt.close()
    return fid_file, fid_txt


def Frame_Length( x, overlap, nwind):
    nx = len(x)
    noverlap = nwind - overlap
    framelen = int((nx - noverlap) / (nwind - noverlap))
    return framelen


def Truelabel2Trueframe( TrueLabel_bin,wsize,wstep ):
    iidx = 0
    Frame_iidx = 0
    Frame_len = Frame_Length(TrueLabel_bin, wstep, wsize)
    Detect = np.zeros([Frame_len, 1])
    while 1 :
        if iidx+wsize <= len(TrueLabel_bin) :
            TrueLabel_frame = TrueLabel_bin[iidx:iidx + wsize - 1]*10
        else:
            TrueLabel_frame = TrueLabel_bin[iidx:]*10

        if (np.sum(TrueLabel_frame) >= wsize / 2) :
            TrueLabel_frame = 1
        else :
            TrueLabel_frame = 0

        if (Frame_iidx >= len(Detect)):
            break

        Detect[Frame_iidx] = TrueLabel_frame
        iidx = iidx + wstep
        Frame_iidx = Frame_iidx + 1
        if (iidx > len(TrueLabel_bin)):
            break

    return Detect


def frame2rawlabel(label, win_len, win_step):

    num_frame = label.shape[0]

    total_len = (num_frame-1) * win_step + win_len
    raw_label = np.zeros((total_len, 1))
    start_indx = 0

    i = 0

    while True:

        if start_indx + win_len > total_len:
            break
        else:
            temp_label = label[i]
            raw_label[
                start_indx+1:start_indx+win_len
            ] = raw_label[start_indx+1:start_indx+win_len] + temp_label
        i += 1

        start_indx = start_indx + win_step

    raw_label = (raw_label >= 1).choose(raw_label, 1)

    return raw_label


def frame2inpt(label, win_len, win_step):
    num_frame = len(label)
    total_len = (num_frame - 1) * win_step + win_len
    raw_label = np.zeros([1, total_len])
    start_indx = 0
    i = 1
    while 1:
        if (start_indx + win_len > total_len):
            break
        temp = label[i] * np.ones([1, win_len])
        raw_label[
            start_indx: start_indx + win_len
        ] = raw_label[start_indx: start_indx + win_len] + temp
        i = i + 1
        start_indx = start_indx + win_step

    raw_label = th_classifier(raw_label, 1)
    return raw_label


def mrcg_extract(audio_dir, save_dir):
    noisy_speech, audio_sr = librosa.load(audio_dir, sr=16000)
    y_label = np.zeros([len(noisy_speech), 1])
    os.makedirs(save_dir, exist_ok=True)
    os.makedirs("{}/Labels".format(save_dir), exist_ok=True)
    name_mrcg = save_dir + '/mrcg'
    name_label = save_dir + '/Labels/label'
    mrcg_mat = np.transpose(mrcg.mrcg_features(noisy_speech, audio_sr))
    winlen = int(np.ceil(audio_sr * 25 * 0.001))
    winstep = int(np.ceil(audio_sr * 10 * 0.001))
    train_mean = np.mean(mrcg_mat, 1)
    train_std = np.std(mrcg_mat, 1)
    framed_label = Truelabel2Trueframe(y_label, winlen, winstep)
    num = 0
    if (len(mrcg_mat) > len(framed_label)) :
        binary_saver(name_mrcg, mrcg_mat[0: len(framed_label),:], num )
        binary_saver(name_label, framed_label, num)
        data_len = len(framed_label)
    else :
        binary_saver(name_mrcg, mrcg_mat, num)
        binary_saver(name_label, framed_label[1: len(mrcg_mat), 1], num )
        data_len = len(mrcg_mat)
    sio.savemat(save_dir+'/normalize_factor',{'train_mean': train_mean, 'train_std': train_std})

    print('MRCG extraction is successifully done.')
    return data_len, winlen, winstep


def vad_func(
    audio_dir,
    mode,
    th,
    output_type,
    is_default,
    off_on_length=20,
    on_off_length=20,
    hang_before=10,
    hang_over=10
):

    os.system('rm -rf result')
    os.system('rm -rf sample_data')
    print('MRCG extraction ...')

    data_len, winlen, winstep = mrcg_extract(audio_dir)

    os.mkdir('./result')

    order = 'python3 ./lib/python/VAD_test.py -m %d -l %d -d %d --data_dir=./sample_data' \
        ' --model_dir=./saved_model --norm_dir=./norm_data' % (mode, data_len, is_default)

    os.system(order)

    pred_result = sio.loadmat('result/pred.mat')
    pp = pred_result['pred']
    result = th_classifier(pp, th)
    result = vad_post(result, off_on_length, on_off_length, hang_before, hang_over)
    if output_type == 1:
        result = frame2rawlabel(result, winlen, winstep)

    return result


def vad_post(post_label, off_on_length=20, on_off_length=20, hang_before=10, hang_over=10):
    '''fill 1 to short valley'''
    offset = False
    onset = False
    for i in range(post_label.shape[0]):

        if i < post_label.shape[0] - 1:
            if post_label[i] == 1 and post_label[i+1] == 0:  # offset detection
                offset = True
                offset_point = i

            if post_label[i] == 0 and post_label[i+1] == 1 and offset:  # offset -> onset detection

                if i - offset_point < off_on_length:
                    post_label[offset_point:i+1, :] = 1  # fill 1 to valley
                    offset = False

    '''remove impulse like detection'''
    post_label = np.concatenate([np.zeros((1, 1)), post_label], axis=0)

    for i in range(post_label.shape[0]):

        if i < post_label.shape[0] - 1:
            if post_label[i] == 0 and post_label[i + 1] == 1:  # onset detection
                onset = True
                onset_point = i

            if post_label[i] == 1 and post_label[i + 1] == 0 and onset:  # onset -> offset detection

                if i - onset_point < on_off_length:
                    post_label[onset_point:i + 1, :] = 0  # fill 0 to hill
                    onset = False

    post_label = post_label[1:]

    '''hang before & over'''

    for i in range(post_label.shape[0]):

        if i < post_label.shape[0] - 1:
            if post_label[i] == 0 and post_label[i + 1] == 1:  # onset detection
                onset = True

                if i - hang_before < 0:
                    post_label[0:i + 1] = 1
                else:
                    post_label[i-hang_before:i + 1] = 1

            if post_label[i] == 1 and post_label[i + 1] == 0 and onset:  # onset -> offset detection

                onset = False
                print(i)
                if i + hang_over > post_label.shape[0]:
                    post_label[i:, :] = 1

                else:
                    post_label[i:i+hang_over, :] = 1

    return post_label