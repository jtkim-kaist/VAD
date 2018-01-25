import tensorflow as tf
import numpy as np
import utils as utils
import re
import data_reader_bDNN_v2 as dr
import os, sys
import time
import subprocess
# import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import scipy.io as sio
from sklearn import metrics
from matplotlib import colors, cm, pyplot as plt
from scipy.stats import binom
from tensorflow.contrib import rnn
from scipy.optimize import brentq
from scipy.interpolate import interp1d

'''program parameter'''

visualization = False
SEED = 1
reset = True  # remove all existed logs and initialize log directories
device = '/gpu:0'
tf.reset_default_graph()
tf.set_random_seed(SEED)
mode = 'test'

# FLAGS = tf.flags.FLAGS
# tf.flags.DEFINE_string('mode', 'test', "mode : train/ test [default : train]")

'''file directory'''

file_dir = "/home/sbie/storage2/VAD_Database/SE_TIMIT_MRCG_0328"
input_dir = file_dir
output_dir = file_dir + "/Labels"
valid_file_dir = "/home/sbie/storage2/VAD_Database/NX_TIMIT_MRCG_small2"
# valid_file_dir = "/home/sbie/storage2/VAD_Database/record_data"

test_file_dir = "/home/sbie/storage2/VAD_Database/NX_TIMIT_MRCG_big"
norm_dir = input_dir
logs_dir = "/home/sbie/github/VAD_bDNN_baseline/logs_proposed_multi"
save_dir = "/home/sbie/storage2/VAD_Database/saved_model/candidate"
initial_logs_dir = "/home/sbie/storage2/VAD_Database/saved_model/my_converted_checkpoint2"
# initial_logs_dir = "/fake_dir"
ckpt_name = '/RF'

if mode == 'test':

    reset = False
    # logs_dir = "/home/sbie/github/VAD_bDNN_baseline/logs_backup0424/logs_proposed930"
    valid_file_dir = "/home/sbie/storage2/VAD_Database/NX_TIMIT_MRCG_small"
    logs_dir = "/home/sbie/storage2/VAD_Database/saved_model/my_converted_checkpoint"
    # test_file_dir = valid_file_dir

if reset:
    print('log directory was initialized')
    os.popen('rm -rf ' + logs_dir + '/*')
    os.popen('mkdir ' + logs_dir + '/train')
    os.popen('mkdir ' + logs_dir + '/valid')

summary_list = ["cost", "accuracy_SNR_-5", "accuracy_SNR_0", "accuracy_SNR_5", "accuracy_SNR_10",
                "accuracy_across_all_SNRs"]

'''in-output parameter'''

w = 19  # w default = 19
u = 9  # u default = 9
assert (w-1) % u == 0, "w-1 must be divisible by u"
num_features = 768  # for MRCG feature
bdnn_winlen = (((w-1) / u) * 2) + 3
bdnn_inputsize = int(bdnn_winlen * num_features)
bdnn_outputsize = int(bdnn_winlen)
model_config = {"w": w, "u": u}

'''training parameter'''

beta = 1  # softmax parameter
SMALL_NUM = 1e-5
eval_num_batches = 2e5
max_epoch = int(5000)
dropout_rate = 0.5
decay = 0.9  # batch normalization decay factor
eval_th = th = 0.5  # 0.5
batch_size = int(4096)  # default = 4096 * 2
valid_batch_size = batch_size
clip_th = 11  # default : 0.90491669
# initLr = 0.000605  # default : 0.000970598, 0.000605
initLr = 0.000605  # default : 0.000970598, 0.000605

lrDecayRate = .95
lrDecayFreq = 200  # default : 200
val_start_step = 100
val_freq = 1
data_len = None
eval_type = 0
'''Model parameter'''

glimpse_hidden = 128
bp_hidden = 128
glimpse_out = bp_out = 128
nGlimpses = 7  # 7
lstm_cell_size = 128
action_hidden_1 = 256  # default : 256
action_hidden_2 = 256  # default : 256

'''attention visualization'''
attention = []


def train_config(c_train_dir, c_valid_dir, c_logs_dir, c_batch_size_eval, c_max_epoch, c_mode):

    global file_dir
    global input_dir
    global output_dir
    global valid_file_dir
    global norm_dir
    global initial_logs_dir
    global logs_dir
    global ckpt_name
    global batch_size
    global valid_batch_size
    global mode
    global max_epoch

    file_dir = c_train_dir
    valid_file_dir = c_valid_dir
    input_dir = file_dir
    output_dir = file_dir + "/Labels"

    norm_dir = file_dir
    initial_logs_dir = logs_dir = c_logs_dir
    # batch_size = valid_batch_size = c_batch_size_eval + 2 * w
    batch_size = valid_batch_size = c_batch_size_eval

    max_epoch = c_max_epoch
    mode = c_mode


def test_config(c_test_dir, c_norm_dir, c_initial_logs_dir, c_batch_size_eval, c_data_len):
    global test_file_dir
    global norm_dir
    global initial_logs_dir
    global ckpt_name
    global valid_batch_size
    global data_len
    global batch_size
    test_file_dir = c_test_dir
    norm_dir = c_norm_dir
    initial_logs_dir = c_initial_logs_dir + '/backup_ckpt'
    print(initial_logs_dir)
    batch_size = valid_batch_size = c_batch_size_eval
    data_len = c_data_len


def config(c_initLr=initLr, c_clip_threshold=clip_th, c_device=device, c_max_epoch=max_epoch):

    global initLr
    global clip_th
    global device
    global max_epoch

    initLr = c_initLr
    clip_th = c_clip_threshold
    device = c_device
    max_epoch = c_max_epoch


def smooth_softmax(x):

    return tf.sigmoid(x) / tf.expand_dims(tf.reduce_sum(tf.sigmoid(x), axis=1), axis=1)


def softmax(x, b):

    return tf.exp(b*x) / tf.expand_dims(tf.reduce_sum(tf.exp(b*x), axis=1), axis=1)


def affine_transform(x, output_dim, seed=0, name=None):
    """
    affine transformation Wx+b
    assumes x.shape = (batch_size, num_features)
    """
    initializer = tf.truncated_normal_initializer(stddev=0.02, seed=seed)

    # weights = tf.get_variable(name + "_w", [x.get_shape()[1], output_dim],
    #                           initializer=tf.contrib.layers.xavier_initializer(seed=seed))
    weights = tf.get_variable(name + "_w", [x.get_shape()[1], output_dim],
                              initializer=initializer)
    b = tf.get_variable(name + "_b", [output_dim], initializer=tf.constant_initializer(0.0))

    return tf.matmul(x, weights) + b


def sw_sensor(inputs, bp):
    global batch_size
    # bp = tf.ones(bp.get_shape().as_list(), dtype=tf.float32) / 7  # for fix the attention
    bp = tf.expand_dims(bp, axis=2)
    bp = tf.tile(bp, (1, 1, num_features))
    # bp = tf.reshape(bp, (inputs.get_shape()[0].value, -1, 1))
    bp = tf.reshape(bp, (batch_size, -1, 1))

    bp = tf.squeeze(bp)

    print(bp.get_shape().as_list()[0])
    print(inputs.get_shape().as_list()[0])

    sw = bp * inputs

    return sw


def get_glimpse(inputs, bp, reuse=None, is_training=True):

    with tf.variable_scope("glimpse_net", reuse=reuse):

        glimpse_input = sw_sensor(inputs, bp)

        act_glimpse_hidden = tf.nn.relu(utils.batch_norm_affine_transform(glimpse_input, glimpse_hidden, decay=decay,
                                                                          name='glimpse_hidden', seed=SEED,
                                                                          is_training=is_training))
        act_bp_hidden = tf.nn.relu(utils.batch_norm_affine_transform(bp, bp_hidden, decay=decay, name='bp_hidden',
                                                                     seed=SEED,
                                                                     is_training=is_training))

        glimpse_feature = tf.nn.relu(utils.batch_norm_affine_transform(act_glimpse_hidden, glimpse_out, decay=decay,
                                                                       name='glimpse_out', seed=SEED,
                                                                       is_training=is_training) +
                                     utils.batch_norm_affine_transform(act_bp_hidden, bp_out, decay=decay,
                                                                       name='bp_out', seed=SEED,
                                                                       is_training=is_training))

    return glimpse_feature


def multinomial_pmf(mean, sample):
    """
    calculate the probability of bernoulli process
    :param mean: mean. shape = (batch_size, num_sbs)
    :param sample: sample. shape = (batch_size, num_sbs)
    :return: p_br: shape = (batch_size, num_sbs)
    """
    p_br = tf.reduce_prod(tf.pow(mean, sample), axis=2)
    return p_br


def bdnn_prediction(batch_size_in, logits, threshold=th):
    bdnn_batch_size = batch_size_in + 2*w
    result = np.zeros((int(bdnn_batch_size), 1))
    indx = np.arange(int(bdnn_batch_size)) + 1
    indx = indx.reshape((int(bdnn_batch_size), 1))
    indx = utils.bdnn_transform(indx, w, u)
    indx = indx[w:(int(bdnn_batch_size)-w), :]
    indx_list = np.arange(w, int(bdnn_batch_size) - w)

    for i in indx_list:
        indx_temp = np.where((indx-1) == i)
        pred = logits[indx_temp]
        pred = np.sum(pred)/pred.shape[0]
        result[i] = pred

    result = np.trim_zeros(result)
    soft_result = np.float32(result)
    result = np.float32(result) >= threshold

    return result.astype(np.float32), soft_result


def summary_generation(eval_file_dir):

    summary_dic = {}

    noise_list = os.listdir(eval_file_dir)
    noise_list = sorted(noise_list)
    summary_dic["summary_ph"] = summary_ph = tf.placeholder(dtype=tf.float32)

    for name in noise_list:

        with tf.variable_scope(name):
            for summary_name in summary_list:
                    summary_dic[name+"_"+summary_name] = tf.summary.scalar(summary_name, summary_ph)

    with tf.variable_scope("Averaged_Results"):

        summary_dic["cost_across_all_noise_types"] = tf.summary.scalar("cost_across_all_noise_types", summary_ph)
        summary_dic["accuracy_across_all_noise_types"]\
            = tf.summary.scalar("accuracy_across_all_noise_types", summary_ph)
        summary_dic["variance_across_all_noise_types"]\
            = tf.summary.scalar("variance_across_all_noise_types", summary_ph)
    return summary_dic


def full_evaluation(m_eval, sess_eval, batch_size_eval, eval_file_dir, summary_writer, summary_dic, itr):

    mean_cost = []
    mean_accuracy = []
    mean_auc = []

    print("-------- Performance for each of noise types --------")

    noise_list = os.listdir(eval_file_dir)
    noise_list = sorted(noise_list)

    summary_ph = summary_dic["summary_ph"]

    for i in range(len(noise_list)):

        noise_name = '/' + noise_list[i]
        eval_input_dir = eval_file_dir + noise_name
        eval_output_dir = eval_file_dir + noise_name + '/Labels'
        eval_data_set = dr.DataReader(eval_input_dir, eval_output_dir, norm_dir, w=w, u=u, name="eval")
        eval_cost, eval_accuracy, eval_list, eval_auc, eval_auc_list = evaluation(m_eval, eval_data_set, sess_eval, batch_size_eval)

        print("--noise type : " + noise_list[i])
        print("cost: %.4f, accuracy across all SNRs: %.4f" % (eval_cost, eval_accuracy*100))

        print('accuracy wrt SNR:')

        print('SNR_-5 : %.4f, SNR_0 : %.4f, SNR_5 : %.4f, SNR_10 : %.4f' % (eval_list[0]*100, eval_list[1]*100,
                                                                            eval_list[2]*100, eval_list[3]*100))
        print('AUC wrt SNR:')
        print('SNR_-5 : %.4f, SNR_0 : %.4f, SNR_5 : %.4f, SNR_10 : %.4f' % (eval_auc_list[0]*100, eval_auc_list[1]*100,
                                                                            eval_auc_list[2]*100, eval_auc_list[3]*100))
        print('')

        eval_summary_list = [eval_cost] + eval_list + [eval_accuracy]

        for j, summary_name in enumerate(summary_list):
            summary_str = sess_eval.run(summary_dic[noise_list[i]+"_"+summary_name], feed_dict={summary_ph: eval_summary_list[j]})
            summary_writer.add_summary(summary_str, itr)

        mean_cost.append(eval_cost)
        mean_accuracy.append(eval_accuracy)
        mean_auc.append(eval_auc)

    mean_cost = np.mean(np.asarray(mean_cost))
    var_accuracy = np.var(np.asarray(mean_accuracy))
    mean_accuracy = np.mean(np.asarray(mean_accuracy))
    mean_auc = np.mean(np.asarray(mean_auc))

    summary_writer.add_summary(sess_eval.run(summary_dic["cost_across_all_noise_types"],
                                             feed_dict={summary_ph: mean_cost}), itr)
    summary_writer.add_summary(sess_eval.run(summary_dic["accuracy_across_all_noise_types"],
                                             feed_dict={summary_ph: mean_accuracy}), itr)
    summary_writer.add_summary(sess_eval.run(summary_dic["variance_across_all_noise_types"],
                                             feed_dict={summary_ph: var_accuracy}), itr)

    print("-------- Performance across all of noise types --------")
    print("cost : %.4f" % mean_cost)
    print("******* averaged accuracy across all noise_types : %.4f *******" % (mean_accuracy*100))
    print("******* averaged auc across all noise_types : %.7f *******" % (mean_auc*100))
    print("******* variance of accuracies across all noise_types : %6.6f *******" % var_accuracy)

    return mean_auc, var_accuracy


def evaluation(m_valid, valid_data_set, sess, eval_batch_size):
    # num_samples = valid_data_set.num_samples
    # num_batches = num_samples / batch_size
    avg_valid_cost = 0.
    avg_valid_accuracy = 0.
    avg_valid_auc = 0.
    avg_sampled_bps = 0.
    itr_sum = 0.

    auc_list = [0 for i in range(valid_data_set._file_len)]
    accuracy_list = [0 for i in range(valid_data_set._file_len)]
    cost_list = [0 for i in range(valid_data_set._file_len)]
    itr_file = 0

    while True:

        valid_inputs, valid_labels = valid_data_set.next_batch(eval_batch_size)

        if valid_data_set.file_change_checker():
            global attention
            auc_list[itr_file] =  avg_valid_auc / itr_sum
            accuracy_list[itr_file] = avg_valid_accuracy / itr_sum
            cost_list[itr_file] = avg_valid_cost / itr_sum
            avg_sampled_bps = avg_sampled_bps / itr_sum
            avg_valid_accuracy = 0.
            avg_valid_cost = 0.
            avg_valid_auc = 0.
            itr_sum = 0
            itr_file += 1
            valid_data_set.file_change_initialize()

            print("%.4f %.4f %.4f %.4f %.4f %.4f %.4f" % (avg_sampled_bps[0], avg_sampled_bps[1], avg_sampled_bps[2],
                                                          avg_sampled_bps[3], avg_sampled_bps[4],
                                                          avg_sampled_bps[5], avg_sampled_bps[6]))
            attention.append(avg_sampled_bps)
            avg_sampled_bps = 0

        if valid_data_set.eof_checker():
            valid_data_set.reader_initialize()
            # print('Valid data reader was initialized!')  # initialize eof flag & num_file & start index
            break

        feed_dict = {m_valid.inputs: valid_inputs, m_valid.labels: valid_labels,
                     m_valid.keep_probability: 1}

        valid_cost, valid_accuracy, valid_soft_result, valid_raw_labels\
            = sess.run([m_valid.cost, m_valid.reward, m_valid.soft_result, m_valid.raw_labels],
                                              feed_dict=feed_dict)

        # auc calculate
        fpr, tpr, thresholds = metrics.roc_curve(valid_raw_labels, valid_soft_result, pos_label=1)
        valid_auc = metrics.auc(fpr, tpr)
        # valid_auc = brentq(lambda x : 1. - x - interp1d(fpr, tpr)(x), 0., 1.)
        # valid_auc = thresholds[np.argmin(abs(tpr-fpr))]
        sampled_bps_tensor = sess.run(m_valid.sampled_bps_tensor, feed_dict=feed_dict)
        avg_sampled_bps += np.mean(sampled_bps_tensor[:, -1, :], axis=0)
        avg_valid_auc += valid_auc
        avg_valid_cost += valid_cost
        avg_valid_accuracy += valid_accuracy
        itr_sum += 1

    total_avg_valid_auc = np.asscalar(np.mean(np.asarray(auc_list)))
    total_avg_valid_cost = np.asscalar(np.mean(np.asarray(cost_list)))
    total_avg_valid_accuracy = np.asscalar(np.mean(np.asarray(accuracy_list)))

    return total_avg_valid_cost, total_avg_valid_accuracy, accuracy_list, total_avg_valid_auc, auc_list


def dense_to_one_hot(labels_dense, num_classes=2):

    """Convert class labels from scalars to one-hot vectors."""
    # copied from TensorFlow tutorial
    num_labels = labels_dense.shape[0]
    index_offset = np.arange(num_labels) * num_classes
    labels_one_hot = np.zeros((num_labels, num_classes))
    labels_one_hot.flat[(index_offset + labels_dense.ravel()).astype(int)] = 1
    return labels_one_hot.astype(np.float32)


class Model(object):

    def __init__(self, batch_size, reuse=None, is_training=True):
        self.cell_outputs = []
        self.batch_size = batch_size
        self.keep_probability = tf.placeholder(tf.float32, name="keep_probabilty")
        self.inputs = tf.placeholder(tf.float32, shape=[None, bdnn_inputsize],
                                              name="inputs")
        self.labels = tf.placeholder(tf.float32, shape=[None, bdnn_outputsize], name="labels")
        self.is_training = is_training
        self.mean_bps = []
        self.sampled_bps = []
        self.baselines = []
        self.global_step = tf.Variable(0, trainable=False)
        self.lr = tf.train.exponential_decay(initLr, self.global_step, lrDecayFreq, lrDecayRate, staircase=True)
        self.raw_reward = 0
        # set inference graph
        cell_outputs = self.inference(reuse)  # (batch_size, bdnn_outputsize)
        # set objective function

        self.cost, self.reward, self.train_op, self.avg_b, self.rminusb, self.sampled_bps_tensor, self.p_bps,\
            self.print_lr, self.soft_result, self.raw_labels = self.calc_reward(cell_outputs)

    def inference(self, reuse=None):

        # initialization
        raw_inputs = self.inputs
        batch_size = self.batch_size
        keep_prob = self.keep_probability
        is_training = self.is_training

        tf.set_random_seed(SEED)  # initialize the random seed at graph level

        lstm_cell = rnn.LayerNormBasicLSTMCell(lstm_cell_size, dropout_keep_prob=keep_prob, reuse=reuse,
                                               dropout_prob_seed=SEED)

        initial_state = lstm_cell.zero_state(batch_size, tf.float32)

        init_sw = tf.ones([batch_size, int(bdnn_winlen)]) * 0  # start sign

        self.mean_bps.append(init_sw)

        init_sw = tf.cast(tf.greater(init_sw, 0.4), tf.float32)
        self.sampled_bps.append(init_sw)

        reuse_recurrent = None

        init_glimpse = self.get_glimpse(raw_inputs, init_sw, reuse=reuse_recurrent)  # (batch_size, glimpse_out)

        inputs = [0] * nGlimpses
        outputs = [0] * nGlimpses
        glimpse = init_glimpse

        for time_step in range(nGlimpses):

            if time_step == 0:
                with tf.variable_scope("core_network", reuse=reuse_recurrent):
                    (cell_output, cell_state) = lstm_cell(glimpse, initial_state)
                    self.cell_outputs.append(initial_state)
            else:
                reuse_recurrent = True
                with tf.variable_scope("core_network", reuse=reuse_recurrent):
                    (cell_output, cell_state) = lstm_cell(glimpse, cell_state)

            inputs[time_step] = glimpse
            outputs[time_step] = cell_output

            if time_step != nGlimpses - 1:  # not final time_step

                glimpse = self.get_next_input(cell_output, reuse=reuse_recurrent)

            else:  # final time_step
                with tf.variable_scope("baseline", reuse=reuse_recurrent):

                    baseline = tf.sigmoid(affine_transform(((cell_output)), 1, name='baseline'))

                    self.baselines.append(baseline)

        return outputs

    def sw_sensor(self, inputs, bp):

        # bp = tf.ones(bp.get_shape().as_list(), dtype=tf.float32) / 7  # for fix the attention
        bp = tf.expand_dims(bp, axis=2)
        bp = tf.tile(bp, (1, 1, num_features))
        # bp = tf.reshape(bp, (inputs.get_shape()[0].value, -1, 1))
        bp = tf.reshape(bp, (self.batch_size, -1, 1))

        bp = tf.squeeze(bp)

        sw = bp * inputs

        return sw

    def get_glimpse(self, inputs, bp, reuse=None):

        is_training = self.is_training

        with tf.variable_scope("glimpse_net", reuse=reuse):
            glimpse_input = self.sw_sensor(inputs, bp)

            act_glimpse_hidden = tf.nn.relu(
                utils.batch_norm_affine_transform(glimpse_input, glimpse_hidden, decay=decay,
                                                  name='glimpse_hidden', seed=SEED, is_training=is_training))
            act_bp_hidden = tf.nn.relu(utils.batch_norm_affine_transform(bp, bp_hidden, decay=decay, name='bp_hidden',
                                                                         seed=SEED, is_training=is_training))

            glimpse_feature = tf.nn.relu(utils.batch_norm_affine_transform(act_glimpse_hidden, glimpse_out, decay=decay,
                                                                           name='glimpse_out',
                                                                           seed=SEED, is_training=is_training) +
                                         utils.batch_norm_affine_transform(act_bp_hidden, bp_out, decay=decay,
                                                                           seed=SEED, name='bp_out', is_training=is_training))
        return glimpse_feature

    def get_next_input(self, cell_output, reuse=None):

        raw_inputs = self.inputs
        is_training = self.is_training

        with tf.variable_scope("baseline", reuse=reuse):

            baseline = tf.sigmoid(affine_transform(((cell_output)), 1, name='baseline'))

            self.baselines.append(baseline)

        with tf.variable_scope("selection_network", reuse=reuse):

            mean_bp = smooth_softmax(
                utils.batch_norm_affine_transform(cell_output, int(bdnn_winlen), decay=decay, name='selection',
                                                  is_training=is_training))
            # mean_bp = softmax(
            #     utils.batch_norm_affine_transform(cell_output, int(bdnn_winlen), decay=decay, name='selection',
            #                                       is_training=is_training), beta)

            self.mean_bps.append(mean_bp)

            # rand_seq = tf.random_uniform(mean_bp.get_shape().as_list(), minval=0, maxval=1, seed=SEED)

            if is_training:
                sampled_bp = tf.multinomial(mean_bp, num_samples=1, seed=SEED)
                sampled_bp = utils.onehot_tensor(sampled_bp, bdnn_winlen)
            else:
                sampled_bp = mean_bp

            sampled_bp = tf.stop_gradient(sampled_bp)

            self.sampled_bps.append(sampled_bp)

        return self.get_glimpse(raw_inputs, mean_bp, reuse=True)

    def action_network(self, outputs):
        is_training = self.is_training
        with tf.variable_scope("action_network"):

            h1_out = tf.nn.relu(utils.batch_norm_affine_transform(outputs, action_hidden_1,
                                                                  decay=decay, name='action_hidden_1',
                                                                  seed=SEED, is_training=is_training))
            h1_out = tf.nn.dropout(h1_out, keep_prob=self.keep_probability, seed=SEED)
            h2_out = tf.nn.relu(utils.batch_norm_affine_transform(h1_out, action_hidden_2,
                                                                  decay=decay, name='action_hidden_2', seed=SEED,
                                                                  is_training=is_training))
            h2_out = tf.nn.dropout(h2_out, keep_prob=self.keep_probability, seed=SEED)

        return h2_out

    def bdnn_prediction(self, logits, threshold):

        batch_size_tensor = tf.constant(self.batch_size, dtype=tf.float32)
        th_tenor = tf.constant(threshold, dtype=tf.float32)

        result, soft_result = tf.py_func(bdnn_prediction, [batch_size_tensor, logits, th_tenor], Tout=[tf.float32, tf.float32])

        return result, soft_result

    @staticmethod
    def np_trim_zeros(x):
        return np.trim_zeros(x)

    def calc_reward(self, outputs):

        batch_size = self.batch_size

        # consider the action at the last time step

        outputs = outputs[-1]
        outputs = tf.reshape(outputs, (batch_size, lstm_cell_size))

        # get the baseline

        b = tf.stack(self.baselines)
        b = tf.tile(b, [1, 1, 1])
        b = tf.reshape(tf.transpose(b, [1, 0, 2]), [batch_size, nGlimpses])
        no_grad_b = tf.stop_gradient(b)

        # get the action

        action_out = self.action_network(outputs)
        logits = tf.sigmoid(affine_transform(action_out, int(bdnn_outputsize), seed=SEED, name="softmax"))
        logits = tf.identity(logits, "logits")
        result, soft_result = self.bdnn_prediction(logits, threshold=th)
        # soft_result = tf.identity(soft_result, 'soft_pred')

        # convert list of tensors to one big tensor

        mean_bps = tf.concat(axis=0, values=self.mean_bps)
        mean_bps = tf.reshape(mean_bps, (nGlimpses, self.batch_size, int(bdnn_winlen)))
        mean_bps = tf.transpose(mean_bps, [1, 0, 2])

        sampled_bps = tf.concat(axis=0, values=self.sampled_bps)
        sampled_bps = tf.reshape(sampled_bps, (nGlimpses, self.batch_size, int(bdnn_winlen)))
        sampled_bps = tf.transpose(sampled_bps, [1, 0, 2])

        # reward for all examples in the batch

        raw_indx = int(np.floor(bdnn_outputsize / 2))
        raw_labels = self.labels[:, raw_indx]
        raw_labels = tf.reshape(raw_labels, shape=(-1, 1))
        raw_labels = tf.identity(raw_labels, 'raw_labels')

        R = tf.cast(tf.equal(result, raw_labels), tf.float32)
        soft_R = tf.stop_gradient(tf.cast(tf.abs(tf.subtract(1 - soft_result, raw_labels)), tf.float32))
        soft_R = tf.reshape(soft_R, (batch_size, 1))
        soft_R = tf.tile(soft_R, [1, nGlimpses])

        # R = tf.cast(tf.abs(tf.subtract(1 - soft_result, raw_labels)), tf.float32)
        R = tf.stop_gradient(R)
        R = tf.reshape(R, (batch_size, 1))
        self.raw_reward = R
        R = tf.tile(R, [1, nGlimpses])
        reward = tf.reduce_mean(R)

        # select the window

        p_bps = multinomial_pmf(mean_bps, sampled_bps)
        p_bps = tf.reshape(p_bps, (self.batch_size, nGlimpses))

        # define the cost function
        sv_part = -tf.square(self.labels - logits)
        rf_part = tf.log(p_bps + SMALL_NUM) * (R - no_grad_b)

        # J = sv_part

        J = tf.concat(axis=1, values=[sv_part, rf_part])  # comment for sv only
        J = tf.reduce_sum(J, 1)
        J = J - tf.reduce_mean(tf.square(R - b), 1)  # comment for sv only
        J = tf.reduce_mean(J, 0)

        # cost = -J

        cost = -tf.reduce_mean(J)
        var_list = tf.trainable_variables()
        grads = tf.gradients(cost, var_list)
        grads, _ = tf.clip_by_global_norm(grads, clip_th)
        optimizer = tf.train.AdamOptimizer(self.lr)
        train_op = optimizer.apply_gradients(zip(grads, var_list), global_step=self.global_step)

        return cost, reward, train_op, tf.reduce_mean(b), tf.reduce_mean(R - b), \
               sampled_bps, tf.reduce_mean(p_bps), self.lr, soft_result, raw_labels


def main(prj_dir=None, model=None, mode=None):

    #                               Configuration Part                       #
    if mode is 'train':

        import path_setting as ps

        set_path = ps.PathSetting(prj_dir, model)
        logs_dir = initial_logs_dir = set_path.logs_dir
        input_dir = set_path.input_dir
        output_dir = set_path.output_dir
        norm_dir = set_path.norm_dir
        valid_file_dir = set_path.valid_file_dir

        sys.path.insert(0, prj_dir+'/configure/ACAM')
        import config as cg

        global initLr, dropout_rate, max_epoch, batch_size, valid_batch_size
        initLr = cg.lr
        dropout_rate = cg.dropout_rate
        max_epoch = cg.max_epoch
        batch_size = valid_batch_size = cg.batch_size

        global w, u
        w = cg.w
        u = cg.u

        global bdnn_winlen, bdnn_inputsize, bdnn_outputsize
        bdnn_winlen = (((w-1) / u) * 2) + 3
        bdnn_inputsize = int(bdnn_winlen * num_features)
        bdnn_outputsize = int(bdnn_winlen)

        global glimpse_hidden, bp_hidden, glimpse_out, bp_out, nGlimpses,\
            lstm_cell_size, action_hidden_1, action_hidden_2

        glimpse_hidden = cg.glimpse_hidden
        bp_hidden = cg.bp_hidden
        glimpse_out = bp_out = cg.glimpse_out
        nGlimpses = cg.nGlimpse  # 7
        lstm_cell_size = cg.lstm_cell_size
        action_hidden_1 = cg.action_hidden_1  # default : 256
        action_hidden_2 = cg.action_hidden_2  # default : 256

    #                               Graph Part                                 #

    mean_acc_list = []
    var_acc_list = []

    print('Mode : ' + mode)
    print("Graph initialization...")
    with tf.device(device):
        with tf.variable_scope("model", reuse=None):
            m_train = Model(batch_size=batch_size, reuse=None, is_training=True)
            # m_train(batch_size)
    with tf.device(device):
        with tf.variable_scope("model", reuse=True):
            m_valid = Model(batch_size=valid_batch_size, reuse=True, is_training=False)

    print("Done")

    #                               Summary Part                               #

    print("Setting up summary op...")
    summary_ph = tf.placeholder(dtype=tf.float32)

    with tf.variable_scope("Training_procedure"):

        cost_summary_op = tf.summary.scalar("cost", summary_ph)
        accuracy_summary_op = tf.summary.scalar("accuracy", summary_ph)

    # train_summary_writer = tf.summary.FileWriter(logs_dir + '/train/', max_queue=4)
    # valid_summary_writer = tf.summary.FileWriter(logs_dir + '/valid/', max_queue=4)
    # summary_dic = summary_generation(valid_file_dir)

    print("Done")

    #                               Model Save Part                            #

    print("Setting up Saver...")
    saver = tf.train.Saver()
    ckpt = tf.train.get_checkpoint_state(initial_logs_dir)
    print("Done")

    #                               Session Part                               #

    sess_config = tf.ConfigProto(allow_soft_placement=True, log_device_placement=False)
    sess_config.gpu_options.allow_growth = True
    sess = tf.Session(config=sess_config)

    if mode is 'train':
        train_summary_writer = tf.summary.FileWriter(logs_dir + '/train/', sess.graph, max_queue=2)
        valid_summary_writer = tf.summary.FileWriter(logs_dir + '/valid/', max_queue=2)

    if ckpt and ckpt.model_checkpoint_path:  # model restore
        print("Model restored...")
        print(initial_logs_dir+ckpt_name)
        if mode is 'train':
            saver.restore(sess, ckpt.model_checkpoint_path)
        else:
            saver.restore(sess, initial_logs_dir+ckpt_name)
            saver.save(sess, initial_logs_dir + "/model_ACAM.ckpt", 0)  # model save

        print("Done")

    else:
        sess.run(tf.global_variables_initializer())  # if the checkpoint doesn't exist, do initialization

    if mode is 'train':
        train_data_set = dr.DataReader(input_dir, output_dir, norm_dir, w=w, u=u,
                                       name="train")  # training data reader initialization
    if mode is 'train':

        for itr in range(max_epoch):

            start_time = time.time()

            train_inputs, train_labels = train_data_set.next_batch(batch_size)

            feed_dict = {m_train.inputs: train_inputs, m_train.labels: train_labels,
                         m_train.keep_probability: dropout_rate}

            sess.run(m_train.train_op, feed_dict=feed_dict)

            if itr % 10 == 0 and itr >= 0:

                train_cost, train_reward, train_avg_b, train_rminusb, train_p_bps, train_lr \
                    = sess.run([m_train.cost, m_train.reward, m_train.avg_b, m_train.rminusb, m_train.p_bps,
                                m_train.print_lr]
                               , feed_dict=feed_dict)

                duration = time.time() - start_time
                print("Step: %d, cost: %.4f, accuracy: %4.4f, b: %4.4f, R-b: %4.4f, p_bps: %4.4f, lr: %7.6f (%.3f sec)"
                      % (itr, train_cost, train_reward, train_avg_b, train_rminusb, train_p_bps, train_lr, duration))

                train_cost_summary_str = sess.run(cost_summary_op, feed_dict={summary_ph: train_cost})
                train_accuracy_summary_str = sess.run(accuracy_summary_op, feed_dict={summary_ph: train_reward})
                train_summary_writer.add_summary(train_cost_summary_str, itr)  # write the train phase summary to event files
                train_summary_writer.add_summary(train_accuracy_summary_str, itr)

            # if train_data_set.eof_checker():

            # if itr % val_freq == 0 and itr >= val_start_step:
            if itr % 50 == 0 and itr > 0:
                saver.save(sess, logs_dir + "/model.ckpt", itr)  # model save
                print('validation start!')
                valid_accuracy, valid_cost = \
                    utils.do_validation(m_valid, sess, valid_file_dir, norm_dir,
                                        type='ACAM')

                print("valid_cost: %.4f, valid_accuracy=%4.4f" % (valid_cost, valid_accuracy * 100))
                valid_cost_summary_str = sess.run(cost_summary_op, feed_dict={summary_ph: valid_cost})
                valid_accuracy_summary_str = sess.run(accuracy_summary_op, feed_dict={summary_ph: valid_accuracy})
                valid_summary_writer.add_summary(valid_cost_summary_str, itr)  # write the train phase summary to event files
                valid_summary_writer.add_summary(valid_accuracy_summary_str, itr)

                # mean_accuracy, var_accuracy = full_evaluation(m_valid, sess, valid_batch_size, valid_file_dir, valid_summary_writer, summary_dic, itr)
                # if mean_accuracy >= 0.991:
                #
                #     print('model was saved!')
                #     model_name = '/model' + str(int(mean_accuracy * 1e4)) + 'and'\
                #                  + str(int(var_accuracy * 1e5)) + '.ckpt'
                #     saver.save(sess, save_dir + model_name, itr)
                # mean_acc_list.append(mean_accuracy)
                # var_acc_list.append(var_accuracy)

                # train_data_set.initialize()

    elif mode == 'test':

        final_softout, final_label = utils.vad_test(m_valid, sess, valid_batch_size, test_file_dir, norm_dir, data_len,
                                                    eval_type)


        # if data_len is None:
        #     return final_softout, final_label
        # else:
        #     final_softout = final_softout[0:data_len, :]
        #     final_label = final_label[0:data_len, :]

        # fpr, tpr, thresholds = metrics.roc_curve(final_label, final_softout, pos_label=1)
        # eval_auc = metrics.auc(fpr, tpr)
        # print(eval_auc)

        # full_evaluation(m_valid, sess, valid_batch_size, test_file_dir, valid_summary_writer, summary_dic, 0)
        # if visualization:
        #     global attention
        #     attention = np.asarray(attention)
        #     sio.savemat('attention.mat', {'attention' : attention})
        #     subprocess.call(['./visualize.sh'])
        if data_len is None:
            return final_softout, final_label
        else:
            return final_softout[0:data_len, :], final_label[0:data_len, :]

if __name__ == "__main__":
    tf.app.run()


