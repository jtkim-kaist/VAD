import tensorflow as tf
import numpy as np
import utils as utils
import re
import data_reader_bDNN as dr
from tensorflow.contrib import rnn

SEED = 1
w = 19  # w default = 19
u = 9  # u default = 9
assert (w-1) % u == 0, "w-1 must be divisible by u"
num_features = 768  # for MRCG feature
bdnn_winlen = (((w-1) / u) * 2) + 3
bdnn_inputsize = int(bdnn_winlen * num_features)
bdnn_outputsize = int(bdnn_winlen)
initLr = 0.000605
# initLr = 0.000605  # default : 0.000970598, 0.000605
lrDecayRate = .95
lrDecayFreq = 200  # default : 200
decay = 0.9  # batch normalization decay factor
rf_threshold = 0.5  # 0.5
SMALL_NUM = 1e-5
clip_th = 11  # default : 0.90491669


glimpse_hidden = 128
bp_hidden = 128
glimpse_out = bp_out = 128
nGlimpses = 7  # 7
lstm_cell_size = 128
action_hidden_1 = 256  # default : 256
action_hidden_2 = 256  # default : 256


class Model(object):

    def __init__(self, batch_size, reuse=None, is_training=True):
        self.cell_outputs = []
        self.batch_size = batch_size
        self.keep_probability = tf.placeholder(tf.float32, name="keep_probabilty")
        self.inputs = tf.placeholder(tf.float32, shape=[batch_size, bdnn_inputsize],
                                              name="inputs")
        self.labels = tf.placeholder(tf.float32, shape=[batch_size, bdnn_outputsize], name="labels")
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
        # lstm_cell = rnn.BasicRNNCell(lstm_cell_size, reuse=reuse)

        initial_state = lstm_cell.zero_state(batch_size, tf.float32)

        init_sw = tf.ones([batch_size, int(bdnn_winlen)]) * 0  # start sign

        self.mean_bps.append(init_sw)

        init_sw = tf.cast(tf.greater(init_sw, 0.4), tf.float32)
        self.sampled_bps.append(init_sw)

        reuse_recurrent = False

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

                    baseline = tf.sigmoid(utils.affine_transform(((cell_output)), 1, name='baseline'))

                    self.baselines.append(baseline)

        return outputs

    def get_glimpse(self, inputs, bp, reuse=None):

        is_training = self.is_training

        with tf.variable_scope("glimpse_net", reuse=reuse):
            glimpse_input = sw_sensor(inputs, bp)

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

            baseline = tf.sigmoid(utils.affine_transform((((cell_output))), 1, name='baseline'))

            self.baselines.append(baseline)

        with tf.variable_scope("selection_network", reuse=reuse):

            mean_bp = smooth_softmax(
                utils.batch_norm_affine_transform((cell_output), int(bdnn_winlen), decay=decay, name='selection',
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

        batch_size_tensor = tf.constant(self.batch_size+2*w, dtype=tf.float32)
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
        logits = tf.sigmoid(utils.affine_transform(action_out, int(bdnn_outputsize), seed=SEED, name="softmax"))
        result, soft_result = self.bdnn_prediction(logits, threshold=rf_threshold)

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
        # rf_part = tf.log(p_bps + SMALL_NUM) * (soft_R - no_grad_b)
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


def smooth_softmax(x):

    return tf.sigmoid(x) / tf.expand_dims(tf.reduce_sum(tf.sigmoid(x), axis=1), axis=1)


def softmax(x, b):

    return tf.exp(b*x) / tf.expand_dims(tf.reduce_sum(tf.exp(b*x), axis=1), axis=1)


def sw_sensor(inputs, bp):

    bp = tf.expand_dims(bp, axis=2)
    bp = tf.tile(bp, (1, 1, num_features))
    bp = tf.reshape(bp, (inputs.get_shape()[0].value, -1, 1))
    bp = tf.squeeze(bp)

    sw = bp * inputs

    return sw


def multinomial_pmf(mean, sample):
    """
    calculate the probability of bernoulli process
    :param mean: mean. shape = (batch_size, num_sbs)
    :param sample: sample. shape = (batch_size, num_sbs)
    :return: p_br: shape = (batch_size, num_sbs)
    """
    p_br = tf.reduce_prod(tf.pow(mean, sample), axis=2)
    return p_br


def bdnn_prediction(bdnn_batch_size, logits, threshold=0.5):

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

def dense_to_one_hot(labels_dense, num_classes=2):

    """Convert class labels from scalars to one-hot vectors."""
    # copied from TensorFlow tutorial
    num_labels = labels_dense.shape[0]
    index_offset = np.arange(num_labels) * num_classes
    labels_one_hot = np.zeros((num_labels, num_classes))
    labels_one_hot.flat[(index_offset + labels_dense.ravel()).astype(int)] = 1
    return labels_one_hot.astype(np.float32)